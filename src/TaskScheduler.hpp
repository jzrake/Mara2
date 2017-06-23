#ifndef TaskScheduler_hpp
#define TaskScheduler_hpp

#include "Timer.hpp"
#include "HDF5.hpp"


class SimulationStatus;
namespace Cow { namespace H5 { class Location; } }



class TaskScheduler
{
public:
	class Recurrence
	{
	public:
		Recurrence (double physTimeInterval, double wallTimeInterval=0.0, int iterationInterval=0);
	private:
		friend class TaskScheduler;
		char isDue (SimulationStatus status, const Cow::Timer& timer) const;
		void update (char reason);
		double physTimeInterval;
		double wallTimeInterval;
		int iterationInterval;
		double nextPhysTime;
		double nextWallTime;
		int nextIteration;
		int repetition;
	};

	class Task
	{
	public:
		virtual void run (SimulationStatus status, int repetition) = 0;
	};

	using TaskFunction = std::function<void (SimulationStatus status, int repetition)>;

    /**
    Add a general task to the scheduler's queue with the given recurrence
    rule.
    */
	void schedule (std::shared_ptr<Task> task, Recurrence recurrence, std::string name="");

    /**
    Add a lambda-function task to the scheduler's queue with the given
    recurrence rule.
    */
	void schedule (TaskFunction task, Recurrence recurrence, std::string name="");

    /**
    Loop over scheduled jobs and run any that are overdue. Jobs are called in
    the order they were scheduled.
    */
	void dispatch (SimulationStatus status);

    /**
    Restore the state of scheduled tasks by reading them from the given HDF5
    location. Only tasks with names are restored.
    */
    void load (Cow::H5::Location& location);

    /**
    Write the state of scheduled tasks to the given HDF5 location. Only tasks
    with names are written.
    */
    void write (Cow::H5::Location& location) const;

private:
	class LambdaTask;
	Cow::Timer timer;
	using RecurrenceTask = std::tuple<std::shared_ptr<Task>, Recurrence, std::string>;
	std::vector<RecurrenceTask> tasks;
};

#endif
