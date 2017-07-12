#include <iomanip>
#include "BlockDecomposition.hpp"
#include "Checkpoint.hpp"
#include "CartesianMeshGeometry.hpp"
#include "FileSystem.hpp"
#include "MeshData.hpp"
#include "TimeSeriesManager.hpp"
#include "TaskScheduler.hpp"

using namespace Cow;

#include "HDF5.hpp"
#include "Timer.hpp"




// ============================================================================
CheckpointWriter::CheckpointWriter()
{
	outputDirectory = "data";
    filenamePrefix = "chkpt";
    format = Format::single;
}

void CheckpointWriter::setFormat (Format formatToUse)
{
    format = formatToUse;
}

void CheckpointWriter::setFormat (std::string formatString)
{
    if (formatString == "single")
    {
        format = Format::single;
        return;
    }
    if (formatString == "multiple")
    {
        format = Format::multiple;
        return;
    }
    throw std::runtime_error ("Checkpoint format string not single or multiple");
}

void CheckpointWriter::setFilenamePrefix (std::string filenamePrefixToUse)
{
    filenamePrefix = filenamePrefixToUse;
}

void CheckpointWriter::setOutputDirectory (std::string outputDirectoryToUse)
{
    outputDirectory = outputDirectoryToUse;
}

void CheckpointWriter::setMeshDecomposition (std::shared_ptr<MeshDecomposition> meshDecompositionToUse)
{
    meshDecomposition = meshDecompositionToUse;
}

void CheckpointWriter::setTimeSeriesManager (std::shared_ptr<TimeSeriesManager> timeSeriesManagerToUse)
{
    timeSeriesManager = timeSeriesManagerToUse;
}

void CheckpointWriter::setTaskScheduler (std::shared_ptr<TaskScheduler> taskSchedulerToUse)
{
    taskScheduler = taskSchedulerToUse;
}

void CheckpointWriter::setUserParameters (Variant::NamedValues userParametersToUse)
{
    userParameters = userParametersToUse;
}

void CheckpointWriter::writeCheckpoint (
	int checkpointNumber,
    SimulationStatus& status,
    ConservationLaw& conservationLaw,
    MeshData& meshData,
    MeshGeometry& meshGeometry,
    Logger& logger) const
{
    auto timer = Timer();
    auto block = dynamic_cast<const BlockDecomposition*> (meshDecomposition.get());
    auto rank = block && (format == Format::multiple) ? block->getCommunicator().rank() : -1;
    auto filename = FileSystem::makeFilename (outputDirectory,
        filenamePrefix, ".h5", checkpointNumber, rank);

    if (meshDecomposition != nullptr && block == nullptr)
    {
    	throw std::runtime_error ("A MeshDecomposition object was given, "
            "but was not a BlockDecomposition");
    }

    auto writeHeaderPart = [&] ()
    {
        // Create data directory, HDF5 checkpoint file, and groups
        // --------------------------------------------------------------------
        FileSystem::ensureDirectoryExists (outputDirectory);
        auto file            = H5::File (filename, "w");
        auto statusGroup     = file.createGroup ("status");
        auto schedulerGroup  = file.createGroup ("scheduler");
        auto timeSeriesGroup = file.createGroup ("time_series");
        auto userGroup       = file.createGroup ("user");
        auto primitiveGroup  = file.createGroup ("primitive");
        auto diagnosticGroup = file.createGroup ("diagnostic");
        auto meshGroup       = file.createGroup ("mesh");


        // Write misc data like date and run script
        // --------------------------------------------------------------------
        auto t = std::time (nullptr);
        auto tm = *std::localtime (&t);
        auto oss = std::ostringstream();
        oss << std::put_time (&tm, "%m-%d-%Y %H:%M:%S");

        file.writeString ("run_name", "Mara");
        file.writeString ("date", oss.str());


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        auto globalShape = Array::vectorFromShape (block && format == Format::single ? block->getGlobalShape() : meshGeometry.cellsShape());
        auto localsShape = Array::vectorFromShape (meshGeometry.cellsShape());
        auto dtype = H5::DataType::nativeDouble();
        auto plist = H5::PropertyList::DataSetCreate().setChunk (localsShape);


        for (int q = 0; q < conservationLaw.getNumConserved(); ++q)
        {
            auto field = conservationLaw.getPrimitiveName(q);
            primitiveGroup.createDataSet (field, globalShape, dtype, plist);
        }

        for (int q = 0; q < meshData.getNumDiagnostics(); ++q)
        {
            auto field = meshData.getDiagnosticName(q);
            diagnosticGroup.createDataSet (field, globalShape, dtype, plist);
        }


        // Write the simulation status data into the checkpoint
        // --------------------------------------------------------------------
        for (auto member : status.pack())
        {
            statusGroup.writeVariant (member.first, member.second);
        }


        // Write the user parameters
        // --------------------------------------------------------------------
        for (auto member : userParameters)
        {
            userGroup.writeVariant (member.first, member.second);
        }


        // Write the mesh information (needs to be generalized to non-rectilinear meshes)
        // --------------------------------------------------------------------
        meshGroup.writeString ("type", "cartesian");

        auto pointGroup = meshGroup.createGroup ("points");
        auto geometry = dynamic_cast<const CartesianMeshGeometry*>
        (block ? block->getGlobalGeometry().get() : &meshGeometry);

        pointGroup.writeArray ("x", geometry->getPointCoordinates (0));
        pointGroup.writeArray ("y", geometry->getPointCoordinates (1));
        pointGroup.writeArray ("z", geometry->getPointCoordinates (2));

        if (block && format == Format::multiple)
        {
            auto S = Array::vectorFromShape (block->getPatchRegion().lower);
            auto T = Array::vectorFromShape (block->getGlobalShape());

            meshGroup.writeVectorInt ("patch_lower", S);
            meshGroup.writeVectorInt ("global_shape", T);
        }

        // Write the time series data and scheduler state
        // --------------------------------------------------------------------
        if (timeSeriesManager) timeSeriesManager->write (timeSeriesGroup);
        if (taskScheduler) taskScheduler->write (schedulerGroup);
    };


    auto writeDataPart = [&] (int rank)
    {
        // Open the file and groups
        // --------------------------------------------------------------------
        auto file            = H5::File (filename, "a");
        auto primitiveGroup  = file.getGroup ("primitive");
        auto diagnosticGroup = file.getGroup ("diagnostic");


        // Create data sets for the primitive and diagnostic data
        // --------------------------------------------------------------------
        auto targetRegion = block && format == Format::single ? block->getPatchRegion() : Region();

        for (int q = 0; q < conservationLaw.getNumConserved(); ++q)
        {
            auto field = conservationLaw.getPrimitiveName(q);
            auto dset = primitiveGroup.getDataSet (field);
            dset[targetRegion] = meshData.getPrimitive(q);
        }

        for (int q = 0; q < meshData.getNumDiagnostics(); ++q)
        {
            auto field = meshData.getDiagnosticName(q);
            auto dset = diagnosticGroup.getDataSet (field);
            dset[targetRegion] = meshData.getDiagnostic(q);
        }
    };



	logger.log ("Checkpoint") << "writing " << filename << std::endl;

    if (block && format == Format::single)
    {
    	auto comm = block->getCommunicator();
	    comm.onMasterOnly (writeHeaderPart);
	    comm.inSequence (writeDataPart);
	}
	else
	{
		writeHeaderPart();
		writeDataPart(0);
	}
    logger.log ("Checkpoint") << "required " << timer.ageInSeconds() << std::endl;
}

void CheckpointWriter::readCheckpoint (
	std::string filename,
    SimulationStatus& status,
    ConservationLaw& conservationLaw,
    MeshData& meshData,
    Logger& logger) const
{
    logger.log ("Checkpoint") << "reading checkpoint file " << filename << std::endl;

    auto block = dynamic_cast<const BlockDecomposition*> (meshDecomposition.get());

    if (meshDecomposition != nullptr && block == nullptr)
    {
        throw std::runtime_error ("A MeshDecomposition object was given, "
            "but was not a BlockDecomposition");
    }

    auto file            = H5::File (filename, "r");
    auto statusGroup     = file.getGroup ("status");
    auto schedulerGroup  = file.getGroup ("scheduler");
    auto primitiveGroup  = file.getGroup ("primitive");
    auto timeSeriesGroup = file.getGroup ("time_series");

    auto sourceRegion = block ? block->getPatchRegion() : Region();
    auto primitive = primitiveGroup.readArrays (conservationLaw.getPrimitiveNames(), 3, sourceRegion);

    meshData.assignPrimitive (primitive);
    status.update (statusGroup.readNamedValues());
    
    if (timeSeriesManager) timeSeriesManager->load (timeSeriesGroup);
    if (taskScheduler) taskScheduler->load (schedulerGroup);
}




// ============================================================================
int CheckpointStitcherProgram::run (int argc, const char* argv[])
{
    auto user = Variant::NamedValues();
    user["outfile"]  = "chkpt.h5";
    Variant::updateFromCommandLine (user, argc, argv);

    auto targetCheckpointName = user["outfile"];
    auto chkpt           = H5::File (targetCheckpointName, "w"); // output file name
    auto chunk           = H5::File (argv[1], "r"); // representative chunk file
    auto primitiveGroup  = chkpt.createGroup ("primitive");
    auto diagnosticGroup = chkpt.createGroup ("diagnostic");
    auto meshGroup       = chkpt.createGroup ("mesh");
    auto globalShape     = chunk.getGroup ("mesh").readVectorInt ("global_shape");
    auto dtype           = H5::DataType::nativeDouble();



    chunk.iterate ([&] (std::string name)
    {
        if (name != "primitive" && name != "diagnostic" && name != "mesh")
        {
            chunk.copy (name, chkpt);
        }
    });

    chunk.getGroup ("mesh").iterate ([&] (std::string name)
    {
        // These data sets are specific to checkpoint files that are only
        // chunks.
        if (name != "global_shape" && name != "patch_lower")
        {
            chunk.getGroup("mesh").copy (name, meshGroup);
        }
    });

    chunk.getGroup ("primitive").iterate ([&] (std::string field)
    {
        primitiveGroup.createDataSet (field, globalShape, H5::DataType::nativeDouble());
    });

    chunk.getGroup ("diagnostic").iterate ([&] (std::string field)
    {
        diagnosticGroup.createDataSet (field, globalShape, H5::DataType::nativeDouble());
    });

    for (int n = 1; n < argc; ++n)
    {
        if (std::string (argv[n]).find ("=") != std::string::npos)
        {
            continue;
        }

        std::cout << "copying " << argv[n] << " -> " << targetCheckpointName << std::endl;

        auto chunk = H5::File (argv[n], "r");
        auto sourcePrim = chunk.getGroup ("primitive");
        auto targetPrim = chkpt.getGroup ("primitive");
        auto sourceDiag = chunk.getGroup ("diagnostic");
        auto targetDiag = chkpt.getGroup ("diagnostic");

        sourcePrim.iterate ([&] (std::string field)
        {
            auto lower = Array::shapeFromVector (chunk.getGroup ("mesh").readVectorInt ("patch_lower"));
            auto shape = Array::shapeFromVector (sourcePrim.getDataSet (field).getSpace().getShape());
            auto targetRegion = Region();

            for (unsigned int n = 0; n < lower.size(); ++n)
            {
                targetRegion.lower[n] = lower[n];
                targetRegion.upper[n] = lower[n] + shape[n];
            }
            targetPrim.getDataSet (field)[targetRegion] = sourcePrim.readArray (field);
        });

        sourceDiag.iterate ([&] (std::string field)
        {
            auto lower = Array::shapeFromVector (chunk.getGroup ("mesh").readVectorInt ("patch_lower"));
            auto shape = Array::shapeFromVector (sourceDiag.getDataSet (field).getSpace().getShape());
            auto targetRegion = Region();

            for (unsigned int n = 0; n < lower.size(); ++n)
            {
                targetRegion.lower[n] = lower[n];
                targetRegion.upper[n] = lower[n] + shape[n];
            }
            targetDiag.getDataSet (field)[targetRegion] = sourceDiag.readArray (field);
        });
    }

    return 0;
}




// ============================================================================
#include <fstream>
#include "VTK.hpp"

int CheckpointToVtkProgram::run (int argc, const char* argv[])
{
    for (int n = 1; n < argc; ++n)
    {
        doFile (argv[n]);
    }
    return 0;
}

void CheckpointToVtkProgram::doFile (std::string filename) const
{
    using namespace Cow;

    auto h5Filename = filename;
    auto doth5 = filename.find (".h5");
    auto vtkFilename = filename.replace (doth5, 5, ".vtk");

    std::cout << h5Filename << " -> " << vtkFilename << std::endl;

    auto file = H5::File (h5Filename, "r");
    auto points = file.getGroup ("mesh").getGroup ("points");
    auto X = points.readArray ("x");
    auto Y = points.readArray ("y");
    auto Z = points.readArray ("z");
    auto cellsShape = Shape {{X.size() - 1, Y.size() - 1, Z.size() - 1, 1, 1}};
    auto primitive = file.getGroup ("primitive");
    auto diagnostic = file.getGroup ("diagnostic");

    auto vtkStream = std::ofstream (vtkFilename);
    auto vtkMesh = VTK::RectilinearGrid (cellsShape);
    vtkMesh.setTitle (file.readString ("run_name"));
    vtkMesh.setUseBinaryFormat (true);
    vtkMesh.setPointCoordinates (X, 0);
    vtkMesh.setPointCoordinates (Y, 1);
    vtkMesh.setPointCoordinates (Z, 2);

    auto velocityNames = std::vector<std::string> {"velocity1", "velocity2", "velocity3"};
    auto magneticNames = std::vector<std::string> {"magnetic1", "magnetic2", "magnetic3"};
    auto currentJNames = std::vector<std::string> {"current1", "current2", "current3"};

    if (primitive.hasDataSet ("density"))
    {
        vtkMesh.addScalarField ("density", primitive.readArray ("density"));
    }
    if (primitive.hasDataSet ("pressure"))
    {
        vtkMesh.addScalarField ("pressure", primitive.readArray ("pressure"));
    }
    if (primitive.hasDataSets (velocityNames))
    {
        vtkMesh.addVectorField ("velocity", primitive.readArrays (velocityNames, 3));
    }
    if (primitive.hasDataSets (magneticNames))
    {
        vtkMesh.addVectorField ("magnetic", primitive.readArrays (magneticNames, 3));
    }
    if (diagnostic.hasDataSets (currentJNames))
    {
        vtkMesh.addVectorField ("current", diagnostic.readArrays (currentJNames, 3));
    }
    if (diagnostic.hasDataSet ("stream_function"))
    {
        vtkMesh.addScalarField ("stream_function", diagnostic.readArray ("stream_function"));
    }

    vtkMesh.write (vtkStream);
}
