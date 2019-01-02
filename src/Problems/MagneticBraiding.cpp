#include <cmath>
#include "MagneticBraiding.hpp"
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
#define SIGN(x) x < 0 ? -1 : 1;

using namespace Cow;




// ============================================================================
int MagneticBraidingProgram::run (int argc, const char* argv[])
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
    user["drive"]   = "orth_shear";
    user["veldr"]   = 0.2;
    user["cool"]    = -1.0;
    user["noise"]   = -1.0;
    user["gamma"]   = 5. / 3;

    auto clineUser = Variant::fromCommandLine (argc, argv);
    auto chkptUser = Variant::NamedValues();

    if (clineUser.find ("restart") != clineUser.end())
    {
        chkptUser = H5::File (clineUser["restart"], "r").getGroup ("user").readNamedValues();
    }

    Variant::update (user, chkptUser);
    Variant::update (user, clineUser);


    auto createBoundaryCondition = [&] ()
    {
        double veldr = user["veldr"];

        auto abc_flceil = [=] (double x, double y, double z)
        {
            const double vx = veldr * std::sin (4 * M_PI * y) * SIGN(z);
            const double vy = veldr * std::cos (4 * M_PI * x) * SIGN(z);
            return std::vector<double> {{vx, vy}};
        };
        auto abc_noceil = [=] (double x, double y, double z)
        {
            const double vx = (z < 0.0) * veldr * std::sin (4 * M_PI * y);
            const double vy = (z < 0.0) * veldr * std::cos (4 * M_PI * x);
            return std::vector<double> {{vx, vy}};
        };
        auto orth_shear = [=] (double x, double y, double z)
        {
            const double a = z > 0.0 ? 1.0 : 0.0;
            const double b = z > 0.0 ? 0.0 : 1.0;
            const double vx = veldr * std::sin (2 * M_PI * y) * a;
            const double vy = veldr * std::cos (2 * M_PI * x) * b;
            return std::vector<double> {{vx, vy}};
        };
        auto single_ft = [=] (double x, double y, double z)
        {
            const double sgnz = SIGN(z);
            const double R = std::sqrt (x * x + y * y);
            const double k = 10.;
            const double v = veldr / k;
            const double vf = v * 3 * k * std::pow (k * R, 2) * std::exp (-std::pow (k * R, 3)); // v-phi;
            const double vx = vf * (-y / R) * sgnz;
            const double vy = vf * ( x / R) * sgnz;
            return std::vector<double> {{vx, vy}};
        };

        auto newBC = new DrivenMHDBoundary;

        if (false) {}
        else if (std::string (user["drive"]) == "abc_flceil") newBC->setBoundaryValueFunction (abc_flceil);
        else if (std::string (user["drive"]) == "abc_noceil") newBC->setBoundaryValueFunction (abc_noceil);
        else if (std::string (user["drive"]) == "orth_shear") newBC->setBoundaryValueFunction (orth_shear);
        else if (std::string (user["drive"]) == "single_ft" ) newBC->setBoundaryValueFunction (single_ft );
        else
        {
            throw std::runtime_error ("unrecognized option for drive");
        }

        return newBC;
    };

    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();
    double timestepSize = 0.0;

    auto bd = std::shared_ptr<BlockDecomposition>();
    auto bc = std::shared_ptr<BoundaryCondition> (createBoundaryCondition());
    auto fs = std::make_shared<MethodOfLinesPlm>();
    auto cs = Shape {{ int (user["N"]), int (user["N"]), int (user["N"]) * int (user["aspect"]) }};
    auto bs = Shape {{ 2, 2, 2 }};
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto cl = std::make_shared<NewtonianMHD>();

    cl->setPressureFloor (1e-4);
    cl->setGammaLawIndex (double (user["gamma"]));
    cl->setCoolingRate (double (user["cool"]));
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
    bc->setMeshGeometry (mg);
    bc->setConservationLaw (cl);
    ct->setMeshSpacing (1.0 / int (user["N"]));

    auto status    = SimulationStatus();
    auto L         = mo->linearCellDimension();
    auto V         = mo->measure (MeshLocation::cell);
    auto P         = md->getPrimitive();
    auto advance   = [&] (double dt) { return ss->advance (*md, 0.0, dt); };
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
        auto Mcell = ct->monopole (Bcell, MeshLocation::cell);
        auto Jcell = ct->current (Bcell, MeshLocation::cell);

        md->allocateDiagnostics ({ "monopole", "current1", "current2", "current3" });
        md->assignDiagnostic (Mcell, 0, MeshData::includeGuard);
        md->assignDiagnostic (Jcell[Region().withRange (3, 0, 1)], 1, MeshData::includeGuard);
        md->assignDiagnostic (Jcell[Region().withRange (3, 1, 2)], 2, MeshData::includeGuard);
        md->assignDiagnostic (Jcell[Region().withRange (3, 2, 3)], 3, MeshData::includeGuard);

        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    };

    auto taskTimeSeries = [&] (SimulationStatus, int rep)
    {
        auto volumeAverageOverPatches = [&] (std::vector<double> vals)
        {
            if (bd)
                return bd->volumeAverageOverPatches (vals);

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
        auto initialData = [&] (double, double, double) -> std::vector<double>
        {
            auto withoutAux = std::vector<double> {1, 0, 0, 0, 1, 0, 0, 1};

            for (int n = 8; n < cl->getNumConserved(); ++n)
            {
                withoutAux.push_back (0.0);
            }

            if (double (user["noise"]) > 0.0)
            {
                withoutAux[0] += double (user["noise"]) * rand() / RAND_MAX;
            }
            return withoutAux;
        };
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    logger->log() << std::endl << user << std::endl;
    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}
