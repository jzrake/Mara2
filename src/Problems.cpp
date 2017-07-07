#include <iomanip>
#include <cmath>
#include "BoundaryConditions.hpp"
#include "CartesianMeshGeometry.hpp"
#include "CellCenteredFieldCT.hpp"
#include "Checkpoint.hpp"
#include "ConservationLaws.hpp"
#include "FieldOperator.hpp"
#include "IntercellFluxSchemes.hpp"
#include "MeshData.hpp"
#include "MeshOperator.hpp"
#include "Problems.hpp"
#include "SolutionSchemes.hpp"
#include "TaskScheduler.hpp"
#include "TimeSeriesManager.hpp"

using namespace Cow;

#include "Timer.hpp"
#include "Logger.cpp"




// ============================================================================
static int maraMainLoop (
    SimulationStatus& status,
    std::function<double ()> timestep,
    std::function<bool ()> condition,
    std::function<void (double)> advance,
    TaskScheduler& scheduler,
    Logger& logger)
{
    auto simulationTimer = Timer();
    auto sessionKzps = 0.0;
    auto sessionIter = 0;

    logger.log ("Mara") << "Beginning main loop" << std::endl;

    while (condition())
    {
        scheduler.dispatch (status);

        auto stepTimer = Timer();
        double dt = timestep();
        advance (dt);

        status.simulationTime += dt;
        status.simulationIter += 1;
        status.wallMinutes = simulationTimer.minutes();

        double kzps = 1e-3 * status.totalCellsInMesh / stepTimer.age();
        logger.log() << "["     << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        logger.log() << "t="    << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        logger.log() << "dt="   << std::setprecision (4) << std::scientific << dt << " ";
        logger.log() << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;

        sessionKzps += kzps;
        sessionIter +=1;
    }

    scheduler.dispatch (status);

    if (sessionIter > 0)
    {
        logger.log ("Mara") << "Completed main loop" << std::endl;
        logger.log ("Mara") << "Mean kzps was " << sessionKzps / sessionIter << std::endl;
    }
    return 0;
}




// ============================================================================
struct Hydro1DTestProgram::Problem
{
    std::string name;
    double finalTime;
    std::shared_ptr<BoundaryCondition> bc;
    InitialDataFunction idf;
    static std::vector<Problem> get();
};

std::vector<Hydro1DTestProgram::Problem> Hydro1DTestProgram::Problem::get()
{
    auto Shocktube1 = [&] (double x, double, double)
    {
        auto S1 = std::vector<double> {1.000, 0.000, 0.0, 0.0, 1.000};
        auto S2 = std::vector<double> {0.125, 0.000, 0.0, 0.0, 0.100};
        return x < 0.5 ? S1 : S2;
    };
    auto Shocktube2 = [&] (double x, double, double)
    {
        auto S1 = std::vector<double> {1.000,-2.000, 0.0, 0.0, 0.400};
        auto S2 = std::vector<double> {1.000, 2.000, 0.0, 0.0, 0.400};
        return x < 0.5 ? S1 : S2;
    };
    auto Shocktube3 = [&] (double x, double, double)
    {
        auto S1 = std::vector<double> {1.0, 0.0, 0.0, 0.0, 1e+3};
        auto S2 = std::vector<double> {1.0, 0.0, 0.0, 0.0, 1e-2};
        return x < 0.5 ? S1 : S2;
    };
    auto Shocktube4 = [&] (double x, double, double)
    {
        auto S1 = std::vector<double> {1.0, 0.0, 0.0, 0.0, 1e-2};
        auto S2 = std::vector<double> {1.0, 0.0, 0.0, 0.0, 1e+2};
        return x < 0.5 ? S1 : S2;
    };
    auto Shocktube5 = [&] (double x, double, double)
    {
        auto S1 = std::vector<double> {5.99924, 19.59750, 0.0, 0.0, 460.894};
        auto S2 = std::vector<double> {5.99924, -6.19633, 0.0, 0.0,  46.095};
        return x < 0.5 ? S1 : S2;
    };
    auto ContactWave = [&] (double x, double, double)
    {
        auto S1 = std::vector<double> {1.0, 0.0, 0.7, 0.2, 1.0};
        auto S2 = std::vector<double> {0.1, 0.0, 0.7, 0.2, 1.0};
        return x < 0.5 ? S1 : S2;
    };
    auto DensityWave = [&] (double x, double, double)
    {
        return std::vector<double> {1.0 + 0.1 * std::sin (4 * M_PI * x), 1.0, 0.0, 0.0, 1.0};
    };

    auto periodic = std::make_shared<PeriodicBoundaryCondition>();
    auto outflow  = std::make_shared<OutflowBoundaryCondition>();
    auto problems = std::vector<Hydro1DTestProgram::Problem>();
    problems.push_back ({ "Shocktube1", 0.100, outflow, Shocktube1 });
    problems.push_back ({ "Shocktube2", 0.100, outflow, Shocktube2 });
    problems.push_back ({ "Shocktube3", 0.005, outflow, Shocktube3 });
    problems.push_back ({ "Shocktube4", 0.01, outflow, Shocktube4 });
    problems.push_back ({ "Shocktube5", 0.01, outflow, Shocktube5 });
    problems.push_back ({ "ContactWave", 0.1, outflow, ContactWave });
    problems.push_back ({ "DensityWave", 1.0, periodic, DensityWave });
    return problems;
}




// ============================================================================
struct Hydro1DTestProgram::Scheme
{
    std::string name;
    std::shared_ptr<SolutionScheme> ss;
    static std::vector<Scheme> get();
};

std::vector<Hydro1DTestProgram::Scheme> Hydro1DTestProgram::Scheme::get()
{
    auto schemes = std::vector<Scheme>();
    auto pcm1 = std::make_shared<MethodOfLinesTVD>();
    auto plm1 = std::make_shared<MethodOfLinesTVD>();
    auto pcm2 = std::make_shared<MethodOfLinesTVD>();
    auto plm2 = std::make_shared<MethodOfLinesTVD>();
    auto pcm3 = std::make_shared<MethodOfLinesTVD>();
    auto plm3 = std::make_shared<MethodOfLinesTVD>();

    pcm1->setIntercellFluxScheme (std::make_shared<MethodOfLines>());
    plm1->setIntercellFluxScheme (std::make_shared<MethodOfLinesPlm>());
    pcm1->setRungeKuttaOrder(1);
    plm1->setRungeKuttaOrder(1);

    pcm2->setIntercellFluxScheme (std::make_shared<MethodOfLines>());
    plm2->setIntercellFluxScheme (std::make_shared<MethodOfLinesPlm>());
    pcm2->setRungeKuttaOrder(2);
    plm2->setRungeKuttaOrder(2);

    pcm3->setIntercellFluxScheme (std::make_shared<MethodOfLines>());
    plm3->setIntercellFluxScheme (std::make_shared<MethodOfLinesPlm>());
    pcm3->setRungeKuttaOrder(3);
    plm3->setRungeKuttaOrder(3);

    schemes.push_back ({ "pcm1", pcm1 });
    schemes.push_back ({ "plm1", plm1 });
    schemes.push_back ({ "pcm2", pcm2 });
    schemes.push_back ({ "plm2", plm2 });
    schemes.push_back ({ "pcm3", pcm3 });
    schemes.push_back ({ "plm3", plm3 });

    return schemes;
}




// ============================================================================
int Hydro1DTestProgram::run (int argc, const char* argv[])
{
    for (const auto& problem : Problem::get())
    {
        for (const auto& scheme : Scheme::get())
        {
            runProblem (problem, scheme);
        }
    }
    return 0;
}

void Hydro1DTestProgram::runProblem (const Problem& problem, const Scheme& scheme)
{
    auto cs = Shape {{128, 1, 1 }}; // cells shape
    auto bs = Shape {{  2, 0, 0 }};
    auto ss = scheme.ss;
    auto bc = problem.bc;
    auto mg = std::make_shared<CartesianMeshGeometry>();
    auto mo = std::make_shared<MeshOperator>();
    auto cl = std::make_shared<NewtonianHydro>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (cs, bs, 5);

    mg->setCellsShape (cs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setFieldOperator (fo);
    ss->setMeshOperator (mo);
    ss->setBoundaryCondition (bc);

    auto status = SimulationStatus();
    auto cflParameter = 0.5;
    auto L = mo->linearCellDimension();
    auto P = md->getPrimitive();

    auto timestep  = [&] ()
    {
        const double dt1 = cflParameter * fo->getCourantTimestep (P, L);
        const double dt2 = problem.finalTime - status.simulationTime;
        return dt1 < dt2 ? dt1 : dt2;
    };
    auto condition = [&] () { return status.simulationTime < problem.finalTime; };
    auto advance   = [&] (double dt) { return ss->advance (*md, dt); };
    auto scheduler = std::make_shared<TaskScheduler>();
    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();

    writer->setFilenamePrefix (problem.name + "-" + scheme.name);
    writer->setMeshDecomposition (nullptr);
    writer->setTimeSeriesManager (nullptr);

    scheduler->schedule ([&] (SimulationStatus, int rep)
    {
        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    }, TaskScheduler::Recurrence (problem.finalTime));

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();

    md->assignPrimitive (mo->generate (problem.idf, MeshLocation::cell));
    md->applyBoundaryCondition (*bc);

    maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);   
}




// ============================================================================
struct Hydro2DTestProgram::Problem
{
    std::string name;
    double finalTime;
    std::shared_ptr<BoundaryCondition> bc;
    InitialDataFunction idf;
    static std::vector<Problem> get();
};

std::vector<Hydro2DTestProgram::Problem> Hydro2DTestProgram::Problem::get()
{
    auto Shocktube2D = [&] (double x, double y, double)
    {
        auto S1 = std::vector<double> {1.000, 0.000, 0.0, 0.0, 1.000};
        auto S2 = std::vector<double> {0.125, 0.000, 0.0, 0.0, 0.100};
        return x + y < 1. ? S1 : S2;
    };
    auto ContactWave2D = [&] (double x, double y, double)
    {
        auto S1 = std::vector<double> {1.0, 0.0, 0.0, 0.0, 1.0};
        auto S2 = std::vector<double> {0.1, 0.0, 0.0, 0.0, 1.0};
        return x + y < 1. ? S1 : S2;
    };
    auto DensityWave2D = [&] (double x, double y, double)
    {
        return std::vector<double> {1.0 + 0.1 * std::sin (4 * M_PI * (x + y)), 1.0, 1.0, 0.0, 1.0};
    };

    auto periodic = std::make_shared<PeriodicBoundaryCondition>();
    auto outflow  = std::make_shared<OutflowBoundaryCondition>();
    auto problems = std::vector<Hydro2DTestProgram::Problem>();
    problems.push_back ({ "Shocktube2D", 0.100, outflow, Shocktube2D });
    problems.push_back ({ "ContactWave2D", 0.1, periodic, ContactWave2D });
    problems.push_back ({ "DensityWave2D", 1.0, periodic, DensityWave2D });
    return problems;
}




// ============================================================================
struct Hydro2DTestProgram::Scheme
{
    std::string name;
    std::shared_ptr<SolutionScheme> ss;
    static std::vector<Scheme> get();
};

std::vector<Hydro2DTestProgram::Scheme> Hydro2DTestProgram::Scheme::get()
{
    auto schemes = std::vector<Scheme>();
    auto pcm2 = std::make_shared<MethodOfLinesTVD>();
    auto plm2 = std::make_shared<MethodOfLinesTVD>();

    pcm2->setIntercellFluxScheme (std::make_shared<MethodOfLines>());
    plm2->setIntercellFluxScheme (std::make_shared<MethodOfLinesPlm>());
    pcm2->setRungeKuttaOrder(2);
    plm2->setRungeKuttaOrder(2);
    schemes.push_back ({ "pcm2", pcm2 });
    schemes.push_back ({ "plm2", plm2 });

    return schemes;
}




// ============================================================================
int Hydro2DTestProgram::run (int argc, const char* argv[])
{
    for (const auto& problem : Problem::get())
    {
        for (const auto& scheme : Scheme::get())
        {
            runProblem (problem, scheme);
        }
    }
    return 0;
}

void Hydro2DTestProgram::runProblem (const Problem& problem, const Scheme& scheme)
{
    auto cs = Shape {{128, 128, 1 }}; // cells shape
    auto bs = Shape {{  2,   2, 0 }};
    auto ss = scheme.ss;
    auto bc = problem.bc;
    auto mg = std::make_shared<CartesianMeshGeometry>();
    auto mo = std::make_shared<MeshOperator>();
    auto cl = std::make_shared<NewtonianHydro>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (cs, bs, 5);

    mg->setCellsShape (cs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setFieldOperator (fo);
    ss->setMeshOperator (mo);
    ss->setBoundaryCondition (bc);

    auto status = SimulationStatus();
    auto cflParameter = 0.5;
    auto L = mo->linearCellDimension();
    auto P = md->getPrimitive();

    auto timestep  = [&] ()
    {
        const double dt1 = cflParameter * fo->getCourantTimestep (P, L);
        const double dt2 = problem.finalTime - status.simulationTime;
        return dt1 < dt2 ? dt1 : dt2;
    };
    auto condition = [&] () { return status.simulationTime < problem.finalTime; };
    auto advance   = [&] (double dt) { return ss->advance (*md, dt); };
    auto scheduler = std::make_shared<TaskScheduler>();
    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();

    writer->setFilenamePrefix (problem.name + "-" + scheme.name);
    writer->setMeshDecomposition (nullptr);
    writer->setTimeSeriesManager (nullptr);

    scheduler->schedule ([&] (SimulationStatus, int rep)
    {
        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    }, TaskScheduler::Recurrence (problem.finalTime));

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();

    md->assignPrimitive (mo->generate (problem.idf, MeshLocation::cell));
    md->applyBoundaryCondition (*bc);

    maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);   
}




// ============================================================================
struct NewtonianMHD2DTestProgram::Problem
{
    std::string name;
    double finalTime;
    std::shared_ptr<BoundaryCondition> bc;
    InitialDataFunction idf;
    InitialDataFunction ivp; // vector potential function
    static std::vector<Problem> get();
};

std::vector<NewtonianMHD2DTestProgram::Problem> NewtonianMHD2DTestProgram::Problem::get()
{
    auto CylindricalBlast = [&] (double x, double y, double)
    {
        x -= 0.5;
        y -= 0.5;
        const double R = std::sqrt (x * x + y * y);
        auto S1 = std::vector<double> {1.000, 0., 0., 0., 1.000, 1.33, 1.33, 0.};
        auto S2 = std::vector<double> {0.125, 0., 0., 0., 0.100, 1.33, 1.33, 0.};
        return R < 0.125 ? S1 : S2;
    };
    auto FieldLoopP = [&] (double x, double y, double)
    {
        x -= 0.5;
        y -= 0.5;
        const double B = 1e-5;
        const double k = 10.0;
        const double R = std::sqrt (x * x + y * y);
        const double bf = B * 3 * k * std::pow (k * R, 2) * std::exp (-std::pow (k * R, 3)); // B-phi
        auto P = std::vector<double>(8);
        P[0] = 1.0;
        P[1] = 1.0;
        P[2] = 1.0;
        P[3] = 0.0;
        P[4] = 1.0;
        P[5] = bf * (-y / R);
        P[6] = bf * ( x / R);
        P[7] = 0.0;
        return P;
    };
    auto FieldLoopA = [&] (double x, double y, double)
    {
        x -= 0.5;
        y -= 0.5;
        const double B = 1e-5;
        const double k = 10.0;
        const double R = std::sqrt (x * x + y * y);
        const double A = B * std::exp (-std::pow (k * R, 3));
        return std::vector<double>  {0.0, 0.0, A };
    };

    const double A = 1.0;
    const double B = 1.0;
    const double C = 0.0;
    const double k = 2 * M_PI;
    auto AbcFieldP = [&] (double x, double y, double z)
    {
        auto bx = C * std::cos (k * z) - B * std::sin (k * y);
        auto by = A * std::cos (k * x) - C * std::sin (k * z);
        auto bz = B * std::cos (k * y) - A * std::sin (k * x);
        return std::vector<double> { 1., 0., 0., 0., 1., bx, by, bz };
    };
    auto AbcFieldA = [&] (double x, double y, double z)
    {
        auto bx = C * std::cos (k * z) - B * std::sin (k * y);
        auto by = A * std::cos (k * x) - C * std::sin (k * z);
        auto bz = B * std::cos (k * y) - A * std::sin (k * x);
        return std::vector<double> { bx / k, by / k, bz / k };
    };

    auto periodic = std::make_shared<PeriodicBoundaryCondition>();
    auto outflow  = std::make_shared<OutflowBoundaryCondition>();
    auto problems = std::vector<NewtonianMHD2DTestProgram::Problem>();
    problems.push_back ({ "CylindricalBlast", 0.1, periodic, CylindricalBlast, nullptr });
    problems.push_back ({ "FieldLoop", 1.0, periodic, FieldLoopP, FieldLoopA });
    problems.push_back ({ "AbcField", 0.1, periodic, AbcFieldP, AbcFieldA });
    return problems;
}




// ============================================================================
struct NewtonianMHD2DTestProgram::Scheme
{
    std::string name;
    std::shared_ptr<SolutionScheme> ss;
    static std::vector<Scheme> get();
};

std::vector<NewtonianMHD2DTestProgram::Scheme> NewtonianMHD2DTestProgram::Scheme::get()
{
    auto schemes = std::vector<Scheme>();
    auto plm2 = std::make_shared<MethodOfLinesTVD>();

    plm2->setIntercellFluxScheme (std::make_shared<MethodOfLinesPlm>());
    plm2->setDisableFieldCT (false);
    plm2->setRungeKuttaOrder(2);
    schemes.push_back ({ "plm2", plm2 });

    return schemes;
}




// ============================================================================
int NewtonianMHD2DTestProgram::run (int argc, const char* argv[])
{
    for (const auto& problem : Problem::get())
    {
        for (const auto& scheme : Scheme::get())
        {
            runProblem (problem, scheme);
        }
    }
    return 0;
}




// ============================================================================
void NewtonianMHD2DTestProgram::runProblem (const Problem& problem, const Scheme& scheme)
{
    auto cs = Shape {{ 64, 64, 1 }}; // cells shape
    auto bs = Shape {{  2,  2, 0 }};
    auto ss = scheme.ss;
    auto bc = problem.bc;
    auto mg = std::make_shared<CartesianMeshGeometry>();
    auto mo = std::make_shared<MeshOperator>();
    auto cl = std::make_shared<NewtonianMHD>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (cs, bs, 8);
    auto ct = std::make_shared<CellCenteredFieldCT>();

    mg->setCellsShape (cs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setFieldOperator (fo);
    ss->setMeshOperator (mo);
    ss->setBoundaryCondition (bc);
    md->setMagneticIndex (cl->getIndexFor (ConservationLaw::VariableType::magnetic));

    auto status = SimulationStatus();
    auto cflParameter = 0.5;
    auto L = mo->linearCellDimension();
    auto P = md->getPrimitive();

    auto timestep  = [&] ()
    {
        const double dt1 = cflParameter * fo->getCourantTimestep (P, L);
        const double dt2 = problem.finalTime - status.simulationTime;
        return dt1 < dt2 ? dt1 : dt2;
    };
    auto condition = [&] () { return status.simulationTime < problem.finalTime; };
    auto advance   = [&] (double dt) { return ss->advance (*md, dt); };
    auto scheduler = std::make_shared<TaskScheduler>();
    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();

    writer->setFilenamePrefix (problem.name + "-" + scheme.name);
    writer->setMeshDecomposition (nullptr);
    writer->setTimeSeriesManager (nullptr);

    scheduler->schedule ([&] (SimulationStatus, int rep)
    {
        auto Bcell = md->getMagneticField (MeshLocation::cell, MeshData::includeGuard);
        auto Mcell = ct->monopole (Bcell, MeshLocation::cell);

        md->allocateDiagnostics ({ "monopole" });
        md->assignDiagnostic (Mcell, 0, MeshData::includeGuard);

        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    }, TaskScheduler::Recurrence (problem.finalTime));

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();

    md->assignPrimitive (mo->generate (problem.idf, MeshLocation::cell));

    if (problem.ivp)
    {
        auto A = mo->generate (problem.ivp, MeshLocation::face, MeshOperator::VectorMode::scalars, bs);
        auto F = ct->vectorPotentialToFluxes (A);
        auto G = ct->generateGodunovFluxes (F, 0);
        auto Bcell = mo->divergence (G);
        md->assignMagneticField (Bcell, MeshLocation::cell, MeshData::includeGuard);
    }
    md->applyBoundaryCondition (*bc);

    maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);   
}




// ============================================================================
#include "BlockDecomposition.hpp"
#define SIGN(x) x < 0 ? -1 : 1;

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
    user["N"]       = 16;
    user["aspect"]  = 1;
    user["drive"]   = "orth_shear";

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
        double vDrive = 0.2;

        auto abc_flceil = [=] (double x, double y, double z)
        {
            const double vx = vDrive * std::sin (4 * M_PI * y) * SIGN(z);
            const double vy = vDrive * std::cos (4 * M_PI * x) * SIGN(z);
            return std::vector<double> {{vx, vy}};
        };
        auto abc_noceil = [=] (double x, double y, double z)
        {
            const double vx = (z < 0.0) * vDrive * std::sin (4 * M_PI * y);
            const double vy = (z < 0.0) * vDrive * std::cos (4 * M_PI * x);
            return std::vector<double> {{vx, vy}};
        };
        auto orth_shear = [=] (double x, double y, double z)
        {
            const double a = z > 0.0 ? 1.0 : 0.0;
            const double b = z > 0.0 ? 0.0 : 1.0;
            const double vx = vDrive * std::sin (2 * M_PI * y) * a;
            const double vy = vDrive * std::cos (2 * M_PI * x) * b;
            return std::vector<double> {{vx, vy}};
        };
        auto single_ft = [=] (double x, double y, double z)
        {
            const double sgnz = SIGN(z);
            const double R = std::sqrt (x * x + y * y);
            const double k = 10.;
            const double v = vDrive / k;
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
    auto cs = Shape {{ int (user["N"]), int (user["N"]), int (user["N"]) * int (user["aspect"]) }};
    auto bs = Shape {{ 2, 2, 2 }};
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);

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
    auto cl = std::make_shared<NewtonianMHD>();
    auto fo = std::make_shared<FieldOperator>();
    auto ct = std::make_shared<CellCenteredFieldCT>();
    auto md = std::make_shared<MeshData> (mg->cellsShape(), bs, 8);

    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setMeshOperator (mo);
    ss->setFieldOperator (fo);
    ss->setBoundaryCondition (bc);
    ss->setRungeKuttaOrder (2);
    ss->setDisableFieldCT (false);
    md->setMagneticIndex (cl->getIndexFor (ConservationLaw::VariableType::magnetic));
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


    if (user["restart"].empty())
    {
        auto initialData = [&] (double, double, double) -> std::vector<double>
        {
            return {1, 0, 0, 0, 1, 0, 0, 1};
        };
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    else
    {
        writer->readCheckpoint (user["restart"], status, *cl, *md, *logger);
        scheduler->skipNext ("checkpoint");
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    logger->log() << std::endl << user << std::endl;
    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}
