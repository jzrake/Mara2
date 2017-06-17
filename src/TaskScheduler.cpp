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
void TaskScheduler::schedule (std::shared_ptr<Task> task, Recurrence recurrence)
{
	tasks.push_back (std::make_pair (recurrence, task));
}

void TaskScheduler::schedule (TaskFunction task, Recurrence recurrence)
{
	schedule (std::make_shared<LambdaTask>(task), recurrence);
}

void TaskScheduler::dispatch (SimulationStatus status)
{
	for (auto& task : tasks)
	{
		if (char reason = task.first.isDue (status, timer))
		{
			task.second->run (status, task.first.repetition);
			task.first.update (reason);
		}
	}
}
