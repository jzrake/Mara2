#include <iomanip>
#include <cmath>
#include "BoundaryConditions.hpp"
#include "CartesianMeshGeometry.hpp"
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
static void maraMainLoop (
    SimulationStatus& status,
    std::function<double ()> timestep,
    std::function<bool ()> condition,
    std::function<void (double)> advance,
    TaskScheduler& scheduler,
    Logger& logger)
{
    auto simulationTimer = Timer();
    auto totKzps = 0.0;
    auto totIter = 0;

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
        logger.log() << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
        logger.log() << "t=" << std::setprecision (4) << std::fixed << status.simulationTime << " ";
        logger.log() << "dt=" << std::setprecision (4) << std::scientific << dt << " ";
        logger.log() << "kzps=" << std::setprecision (2) << std::fixed << kzps << std::endl;

        totKzps += kzps;
        totIter +=1;
    }

    scheduler.dispatch (status);
    logger.log ("Mara") << "Completed main loop" << std::endl;
    logger.log ("Mara") << "Mean kzps was " << totKzps / totIter << std::endl;
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
            // if (problem.name == "DensityWave")
            run (problem, scheme);
        }
    }
    return 0;
}

void Hydro1DTestProgram::run (const Problem& problem, const Scheme& scheme)
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
        return x + y < 0.5 ? S1 : S2;
    };
    auto ContactWave2D = [&] (double x, double y, double)
    {
        auto S1 = std::vector<double> {1.0, 0.0, 0.7, 0.2, 1.0};
        auto S2 = std::vector<double> {0.1, 0.0, 0.7, 0.2, 1.0};
        return x + y < 0.5 ? S1 : S2;
    };
    auto DensityWave2D = [&] (double x, double y, double)
    {
        return std::vector<double> {1.0 + 0.1 * std::sin (4 * M_PI * (x + y)), 1.0, .0, 0.0, 1.0};
    };

    auto periodic = std::make_shared<PeriodicBoundaryCondition>();
    auto outflow  = std::make_shared<OutflowBoundaryCondition>();
    auto problems = std::vector<Hydro2DTestProgram::Problem>();
    problems.push_back ({ "Shocktube2D", 0.100, outflow, Shocktube2D });
    problems.push_back ({ "ContactWave2D", 0.1, outflow, ContactWave2D });
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
            run (problem, scheme);
        }
    }
    return 0;
}

void Hydro2DTestProgram::run (const Problem& problem, const Scheme& scheme)
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
