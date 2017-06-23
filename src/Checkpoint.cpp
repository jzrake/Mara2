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
}

void CheckpointWriter::setOutputDirectory (std::string outputDirectoryToUse)
{
	outputDirectory = outputDirectoryToUse;
}

void CheckpointWriter::setFilenamePrefix (std::string filenamePrefixToUse)
{
    filenamePrefix = filenamePrefixToUse;
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
    auto filename = FileSystem::makeFilename (outputDirectory, filenamePrefix, ".h5", checkpointNumber);
    auto block = dynamic_cast<const BlockDecomposition*> (meshDecomposition.get());

    if (meshDecomposition != nullptr && block == nullptr)
    {
    	throw std::runtime_error ("A MeshDecomposition object was given, but was not a BlockDecomposition");
    }

    auto writeHeaderPart = [&] ()
    {
        // Create data directory, HDF5 checkpoint file, and groups
        // --------------------------------------------------------------------
        FileSystem::ensureDirectoryExists (outputDirectory);
        auto file            = H5::File (filename, "w");
        auto statusGroup     = file.createGroup ("status");
        auto schedulerGroup  = file.createGroup ("scheduler");
        auto userGroup       = file.createGroup ("user");
        auto primitiveGroup  = file.createGroup ("primitive");
        auto diagnosticGroup = file.createGroup ("diagnostic");
        auto meshGroup       = file.createGroup ("mesh");
        auto timeSeriesGroup = file.createGroup ("time_series");


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
        auto globalShape = Array::vectorFromShape (block ? block->getGlobalShape() : meshGeometry.cellsShape());
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

        auto geometry = dynamic_cast<const CartesianMeshGeometry*> (block ? block->getGlobalGeometry().get() : &meshGeometry);
        auto pointGroup = meshGroup.createGroup ("points");
        pointGroup.writeArray ("x", geometry->getPointCoordinates (0));
        pointGroup.writeArray ("y", geometry->getPointCoordinates (1));
        pointGroup.writeArray ("z", geometry->getPointCoordinates (2));


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
        auto targetRegion = block ? block->getPatchRegion() : Region();

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

    if (block)
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
        throw std::runtime_error ("A MeshDecomposition object was given, but was not a BlockDecomposition");
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
