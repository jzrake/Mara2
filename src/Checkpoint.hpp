#ifndef Checkpoint_hpp
#define Checkpoint_hpp

#include "Mara.hpp"

class TaskScheduler;




class CheckpointWriter
{
public:
	CheckpointWriter();

	void setOutputDirectory   (std::string outputDirectoryToUse);
	void setFilenamePrefix    (std::string filenamePrefixToUse);
	void setMeshDecomposition (std::shared_ptr<MeshDecomposition> meshDecompositionToUse);
	void setTimeSeriesManager (std::shared_ptr<TimeSeriesManager> timeSeriesManagerToUse);
    void setTaskScheduler     (std::shared_ptr<TaskScheduler> taskSchedulerToUse);
    void setUserParameters    (Variant::NamedValues userParametersToUse);

	void writeCheckpoint (
		int checkpointNumber,
	    SimulationStatus& status,
	    ConservationLaw& conservationLaw,
	    MeshData& meshData,
	    MeshGeometry& meshGeometry,
	    Logger& logger) const;

    void readCheckpoint (
        std::string filename,
        SimulationStatus& status,
        ConservationLaw& conservationLaw,
        MeshData& meshData,
        Logger& logger) const;

private:
    Variant::NamedValues userParameters;
	std::string filenamePrefix;
	std::string outputDirectory;
	std::shared_ptr<MeshDecomposition> meshDecomposition;
	std::shared_ptr<MeshGeometry>      meshGeometry;
    std::shared_ptr<TimeSeriesManager> timeSeriesManager;
    std::shared_ptr<TaskScheduler>     taskScheduler;
};


#endif
