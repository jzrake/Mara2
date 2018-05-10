#include <cmath>
#include "AlfvenWaveCollision.hpp"
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
int AlfvenWaveCollision::run (int argc, const char* argv[])
{
    auto user = Variant::NamedValues();
    user["outdir"]  = "data";
    user["restart"] = "";
    user["tfinal"]  = 0.1;
    user["serial"]  = false;
    user["cpi"]     = 0.1;
    user["cpf"]     = "single"; // or multiple
    user["tsi"]     = 0.1;
    user["cfl"]     = 0.6;
    user["plm"]     = 1.5;
    user["phi"]     = 0.0; // relative polatization angle (in degrees)
    user["sigma"]   = 1.0; // sigma (B^2 / rho)
    user["A"]       = 0.5; // dB / B
    user["version"] = 1; // IC version: 1 or 2
    user["N"]       = 128;
    user["gamma"]   = 4. / 3;

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
    auto cs = Shape {{ 1, 1, int (user["N"]) }};
    auto bs = Shape {{ 0, 0, 2 }};
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto cl = std::make_shared<RelativisticMHD>();
    // auto cl = std::make_shared<NewtonianMHD>();

    cl->setPressureFloor (1e-4);
    cl->setGammaLawIndex (double (user["gamma"]));
    mg->setCellsShape (cs);
    mg->setLowerUpper ({{-0.5, -0.5, 0.0}}, {{0.5, 0.5, 1.0}});

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
    ss->setRungeKuttaOrder (3);
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
    else if (int (user["version"]) == 1)
    {
        auto initialData = [&] (double, double, double z) -> std::vector<double>
        {
            const double A = double (user["A"]);
            const double dz = 0.1;
            const double d0 = 1.0 / double (user["sigma"]);
            const double p0 = 1.0 / double (user["sigma"]);
            const double Bz = 1.0;
            const double mu = 4. * p0 / d0;
            const double sg = Bz * Bz / d0;
            const double h0 = 1.0 + mu + sg;
            const double va = Bz / std::sqrt (d0 * h0);
            const double phi = M_PI / 180 * double (user["phi"]);

            const double zp = 0.25;
            const double Bp = +std::exp (-std::pow (z - zp, 2) / dz / dz) * A;
            const double vp = -Bp / Bz * va;

            const double zm = 0.75;
            const double Bm = +std::exp (-std::pow (z - zm, 2) / dz / dz) * A;
            const double vm = +Bm / Bz * va;

            double vx = vp + std::cos (phi) * vm;
            double vy = 0  + std::sin (phi) * vm;
            const double Bx = Bp + std::cos (phi) * Bm;
            const double By = 0  + std::sin (phi) * Bm;

            const double g0 = std::sqrt (1.0 + (vx * vx + vy * vy));
            vx /= g0;
            vy /= g0;

            return std::vector<double> {d0, vx, vy, 0, p0, Bx, By, Bz};
        };
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    else if (int (user["version"]) == 2)
    {
        auto initialData = [&] (double, double, double z) -> std::vector<double>
        {
            const double A = double (user["A"]);
            const double dz = 0.1;
            const double d0 = 1.0 / double (user["sigma"]);
            const double p0 = 1.0 / double (user["sigma"]);
            const double Bz = 1.0;
            const double mu = 4. * p0 / d0;
            const double sg = Bz * Bz / d0;
            const double h0 = 1.0 + mu + sg;
            const double va = Bz / std::sqrt (d0 * h0);
            const double phi = M_PI / 180 * double (user["phi"]);

            const double zp = 0.25;
            const double Bp = +std::exp (-std::pow (z - zp, 2) / dz / dz) * A;
            const double vp = -Bp * Bz / (Bz * Bz + Bp * Bp);

            const double zm = 0.75;
            const double Bm = +std::exp (-std::pow (z - zm, 2) / dz / dz) * A;
            const double vm = +Bm * Bz / (Bz * Bz + Bm * Bm);

            double vx = vp + std::cos (phi) * vm;
            double vy = 0  + std::sin (phi) * vm;
            double vz = +Bp * Bp / (Bz * Bz + Bp * Bp) - Bm * Bm / (Bz * Bz + Bm * Bm);
            const double Bx = Bp + std::cos (phi) * Bm;
            const double By = 0  + std::sin (phi) * Bm;

            const double g0 = std::sqrt (1.0 + (vx * vx + vy * vy + vz * vz));
            vx /= g0;
            vy /= g0;
            vz /= g0;

            return std::vector<double> {d0, vx, vy, vz, p0, Bx, By, Bz};
            // return std::vector<double> {d0, vx, 0.0, vz, p0, Bx, 0.0, Bz};
        };
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    logger->log() << std::endl << user << std::endl;
    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}
