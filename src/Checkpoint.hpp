#ifndef Checkpoint_hpp
#define Checkpoint_hpp

#include "Mara.hpp"




class CheckpointWriter
{
public:
	CheckpointWriter();

	void setOutputDirectory (std::string outputDirectoryToUse);
	void setFilenamePrefix (std::string filenamePrefixToUse);
	void setMeshDecomposition (std::shared_ptr<MeshDecomposition> meshDecompositionToUse);
	void setTimeSeriesManager (std::shared_ptr<TimeSeriesManager> timeSeriesManagerToUse);

	void writeCheckpoint (
		int checkpointNumber,
	    SimulationStatus& status,
	    ConservationLaw& conservationLaw,
	    MeshData& meshData,
	    MeshGeometry& meshGeometry,
	    Logger& logger) const;

private:
	std::string filenamePrefix;
	std::string outputDirectory;
	std::shared_ptr<MeshDecomposition> meshDecomposition;
	std::shared_ptr<MeshGeometry> meshGeometry;
	std::shared_ptr<TimeSeriesManager> timeSeriesManager;
};




class BlockDecomposition; // Here temporarily

class CheckpointReader
{
public:
	void readCheckpoint (
		std::string filename,
	    SimulationStatus& status,
        ConservationLaw& conservationLaw,
	    MeshData& meshData,
	    BlockDecomposition& block,
	    TimeSeriesManager& tseries,
	    Logger& logger) const;
};

#endif
