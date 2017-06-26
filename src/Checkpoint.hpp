#ifndef Checkpoint_hpp
#define Checkpoint_hpp

#include "Mara.hpp"

class TaskScheduler;
class BlockDecomposition;



class CheckpointWriter
{
public:
    enum class Format { single, multiple };

	CheckpointWriter();

    void setFormat            (Format formatToUse);
    void setFormat            (std::string formatString);
    void setFilenamePrefix    (std::string filenamePrefixToUse);
    void setOutputDirectory   (std::string outputDirectoryToUse);
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
    Format format;
    Variant::NamedValues userParameters;
    std::string filenamePrefix;
	std::string outputDirectory;
	std::shared_ptr<MeshDecomposition> meshDecomposition;
	std::shared_ptr<MeshGeometry>      meshGeometry;
    std::shared_ptr<TimeSeriesManager> timeSeriesManager;
    std::shared_ptr<TaskScheduler>     taskScheduler;
};



class CheckpointStitcherProgram : public SubProgram
{
public:
    int run (int argc, const char* argv[]) override;
};

#endif
