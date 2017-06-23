#include "Mara.hpp"
#include "TaskScheduler.hpp"




// ============================================================================
TaskScheduler::Recurrence::Recurrence (double physTimeInterval, double wallTimeInterval, int iterationInterval) :
physTimeInterval (physTimeInterval),
wallTimeInterval (wallTimeInterval),
iterationInterval (iterationInterval)
{
	nextPhysTime = 0.0;
	nextWallTime = 0.0;
	nextIteration = 0;
	repetition = 0;
}

char TaskScheduler::Recurrence::isDue (SimulationStatus status, const Cow::Timer& timer) const
{
	if (physTimeInterval > 0.0 && nextPhysTime <= 1e-12 + status.simulationTime) return 'P';
	if (wallTimeInterval > 0.0 && nextWallTime <= timer.minutes()) return 'W';
	if (iterationInterval > 0 && nextIteration <= status.simulationIter) return 'I';
	return '\0';
}

void TaskScheduler::Recurrence::update (char reason)
{
	repetition += 1;

	switch (reason)
	{
		case 'P': nextPhysTime += physTimeInterval; break;
		case 'W': nextWallTime += wallTimeInterval; break;
		case 'I': nextIteration += iterationInterval; break;
	}
}




// ============================================================================
class TaskScheduler::LambdaTask : public TaskScheduler::Task
{
public:
	LambdaTask (TaskFunction taskFunction) : taskFunction (taskFunction) {}
    void run (SimulationStatus status, int rep) override
    {
    	taskFunction (status, rep);
    }
private:
    TaskFunction taskFunction;
};




// ============================================================================
void TaskScheduler::schedule (std::shared_ptr<Task> task, Recurrence recurrence, std::string name)
{
	tasks.push_back (std::make_tuple (task, recurrence, name));
}

void TaskScheduler::schedule (TaskFunction task, Recurrence recurrence, std::string name)
{
	schedule (std::make_shared<LambdaTask>(task), recurrence, name);
}

void TaskScheduler::dispatch (SimulationStatus status)
{
	for (auto& taskDescription : tasks)
	{
        auto& task = std::get<0> (taskDescription);
        auto& recr = std::get<1> (taskDescription);

		if (char reason = recr.isDue (status, timer))
		{
			task->run (status, recr.repetition);
			recr.update (reason);
		}
	}
}

void TaskScheduler::load (Cow::H5::Location& location)
{

}

void TaskScheduler::write (Cow::H5::Location& location) const
{
    for (auto& taskDescription : tasks)
    {
        auto& recr = std::get<1> (taskDescription);
        auto& name = std::get<2> (taskDescription);

        if (! name.empty())
        {
            auto group = location.createGroup (name);
            group.writeDouble ("nextPhysTime", recr.nextPhysTime);
            group.writeDouble ("nextWallTime", recr.nextWallTime);
            group.writeInt ("nextIteration", recr.nextIteration);
            group.writeInt ("repetition", recr.repetition);
        }
    }
}
