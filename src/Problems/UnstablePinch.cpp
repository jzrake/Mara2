#include <cmath>
#include "UnstablePinch.hpp"
#include "../Mara.hpp"
#include "../BoundaryConditions.hpp"
#include "../BlockDecomposition.hpp"
#include "../CartesianMeshGeometry.hpp"
#include "../CellCenteredFieldCT.hpp"
#include "../Checkpoint.hpp"
#include "../ConservationLaws.hpp"
#include "../FieldOperator.hpp"
#include "../IntercellFluxSchemes.hpp"
#include "../MeshData.hpp"
#include "../MeshOperator.hpp"
#include "../Problems.hpp"
#include "../SolutionSchemes.hpp"
#include "../TaskScheduler.hpp"
#include "../TimeSeriesManager.hpp"
#include "MPI.hpp"
#include "HDF5.hpp"
#include "Timer.hpp"
#include "Logger.hpp"

using namespace Cow;


// ============================================================================
int UnstablePinchProgram::run (int argc, const char* argv[])
{
    auto user = Variant::NamedValues();
    user["outdir"]  = "data";
    user["restart"] = "";
    user["tfinal"]  = 16.0;
    user["serial"]  = false;
    user["cpi"]     = 0.25;
    user["cpf"]     = "single"; // or multiple
    user["tsi"]     = 0.1;
    user["cfl"]     = 0.3;
    user["plm"]     = 1.5;
    user["N"]       = 16;
    user["aspect"]  = 1;
    user["L"]       = 0.1;
    user["gamma"]   = 5. / 3;

    auto clineUser = Variant::fromCommandLine (argc, argv);
    auto chkptUser = Variant::NamedValues();

    if (clineUser.find ("restart") != clineUser.end())
    {
        chkptUser = H5::File (clineUser["restart"], "r").getGroup ("user").readNamedValues();
    }
    Variant::update (user, chkptUser);
    Variant::update (user, clineUser);

    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();
    double timestepSize = 0.0;

    auto bd = std::shared_ptr<BlockDecomposition>();
    auto bc = std::shared_ptr<BoundaryCondition> (new PeriodicBoundaryCondition);
    auto fs = std::make_shared<MethodOfLinesPlm>();
    auto cs = Shape {{ int (user["N"]), int (user["N"]), int (user["N"]) * int (user["aspect"]) }};
    auto bs = Shape {{ 2, 2, 2 }};
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    //auto cl = std::make_shared<NewtonianMHD>();
    auto cl = std::make_shared<RelativisticMHD>();

    cl->setPressureFloor (1e-4);
    cl->setGammaLawIndex (double (user["gamma"]));
    cl->setCoolingRate (-1.0);
    mg->setCellsShape (cs);
    mg->setLowerUpper ({{-0.5, -0.5, -0.5 * int (user["aspect"])}}, {{0.5, 0.5, 0.5 * int (user["aspect"])}});

    if (! user["serial"])
    {
        logger->setLogToNullUnless (MpiCommunicator::world().isThisMaster());
        bd = std::make_shared<BlockDecomposition> (mg, *logger);
        mg = bd->decompose();
        bc = bd->createBoundaryCondition (bc);
    }
    tseries->setLogger (logger);
    scheduler->setLogger (logger);

    auto ss = std::make_shared<MethodOfLinesTVD>();
    auto mo = std::make_shared<MeshOperator>();
    auto fo = std::make_shared<FieldOperator>();
    auto ct = std::make_shared<CellCenteredFieldCT>();
    auto md = std::make_shared<MeshData> (mg->cellsShape(), bs, cl->getNumConserved());

    fs->setPlmTheta (double (user["plm"]));
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setMeshOperator (mo);
    ss->setFieldOperator (fo);
    ss->setBoundaryCondition (bc);
    ss->setRungeKuttaOrder (2);
    ss->setDisableFieldCT (false);
    ss->setIntercellFluxScheme (fs);
    md->setMagneticIndex (cl->getIndexFor (ConservationLaw::VariableType::magnetic));
    md->setVelocityIndex (cl->getIndexFor (ConservationLaw::VariableType::velocity));
    bc->setMeshGeometry (mg);
    bc->setConservationLaw (cl);
    ct->setMeshSpacing (1.0 / int (user["N"]));

    auto status    = SimulationStatus();
    auto L         = mo->linearCellDimension();
    auto V         = mo->measure (MeshLocation::cell);
    auto P         = md->getPrimitive();
    auto advance   = [&] (double dt) { return ss->advance (*md, dt); };
    auto condition = [&] () { return status.simulationTime < double (user["tfinal"]); };
    auto timestep  = [&] () { return timestepSize; };

    auto taskRecomputeDt = [&] (SimulationStatus, int rep)
    {
        double localDt = double (user["cfl"]) * fo->getCourantTimestep (P, L);
        timestepSize = MpiCommunicator::world().minimum (localDt);
    };

    auto taskCheckpoint = [&] (SimulationStatus, int rep)
    {
        auto Bcell = md->getMagneticField (MeshLocation::cell, MeshData::includeGuard);
        auto Ecell = md->getElectricField (MeshLocation::cell, MeshData::includeGuard);
        auto Mcell = ct->monopole (Bcell, MeshLocation::cell);
        auto Jcell = ct->current (Bcell, MeshLocation::cell);

        int n = 0;
        md->allocateDiagnostics ({
            "monopole",
            "electric1", "electric2", "electric3" });

        md->assignDiagnostic (Mcell, n++, MeshData::includeGuard);
        md->assignDiagnostic (Ecell[Region().withRange (3, 0, 1)], n++, MeshData::includeGuard);
        md->assignDiagnostic (Ecell[Region().withRange (3, 1, 2)], n++, MeshData::includeGuard);
        md->assignDiagnostic (Ecell[Region().withRange (3, 2, 3)], n++, MeshData::includeGuard);
        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    };

    auto taskTimeSeries = [&] (SimulationStatus, int rep)
    {
        auto volumeAverageOverPatches = [&] (std::vector<double> vals)
        {
            if (bd)
            {
                return bd->volumeAverageOverPatches (vals);
            }

            for (auto& val : vals)
                val /= mg->meshVolume();

            return vals;
        };

        auto volumeIntegrated = fo->volumeIntegratedDiagnostics (P, V);
        auto volumeAveraged = volumeAverageOverPatches (volumeIntegrated);
        auto fieldNames = cl->getDiagnosticNames();
        auto entry = Variant::NamedValues();

        for (unsigned int n = 0; n < fieldNames.size(); ++n)
        {
            entry[fieldNames[n]] = volumeAveraged[n];
        }
        tseries->append (status, entry);
    };


    scheduler->schedule (taskTimeSeries, TaskScheduler::Recurrence (user["tsi"]), "time_series");
    scheduler->schedule (taskCheckpoint, TaskScheduler::Recurrence (user["cpi"]), "checkpoint");
    scheduler->schedule (taskRecomputeDt, TaskScheduler::Recurrence (0.0, 0.0, 16), "compute_dt");

    writer->setMeshDecomposition (bd);
    writer->setTimeSeriesManager (tseries);
    writer->setTaskScheduler     (scheduler);
    writer->setOutputDirectory   (user["outdir"]);
    writer->setFormat            (user["cpf"]);
    writer->setUserParameters    (user);
    writer->setFilenamePrefix    ("chkpt");

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();


    if (! user["restart"].empty())
    {
        writer->readCheckpoint (user["restart"], status, *cl, *md, *logger);
        scheduler->skipNext ("checkpoint");
    }
    else
    {
        auto initialData = [&] (double x, double y, double) -> std::vector<double>
        {
            const double L = user["L"];
            const double R = std::sqrt (x * x + y * y) / L;
            // const double B0 = 5.0; // Values from first round of (non-relativistic) runs
            // const double P0 = 0.1;
            // const double dc = 0.2;
            // const double d0 = 0.1;
            const double B0 = 1.0;
            const double P0 = 1.0;
            const double dc = 2.0;
            const double d0 = 1.0;
            const double Bf = B0 * R * std::exp (1 - R);
            const double Bx = -Bf * y / R;
            const double By = +Bf * x / R;
            const double Pg = P0 + B0 * B0 * std::exp(2) / 4 * (1 - std::exp (-2 * R) * (1 + 2 * R * (1 - R)));
            const double dg = d0 + (dc - d0) / std::pow (std::cosh (2 * R), 2) + (1e-2 * rand()) / RAND_MAX;
            return std::vector<double> {dg, 0, 0, 0, Pg, Bx, By, 0};
        };
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    logger->log() << std::endl << user << std::endl;
    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}