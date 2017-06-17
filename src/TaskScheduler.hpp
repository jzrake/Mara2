#ifndef TaskScheduler_hpp
#define TaskScheduler_hpp

#include "Timer.hpp"

class SimulationStatus;




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

	void schedule (std::shared_ptr<Task> task, Recurrence recurrence);
	void dispatch (SimulationStatus status);

private:
	Cow::Timer timer;
	using RecurrenceTask = std::pair<Recurrence, std::shared_ptr<Task>>;
	std::vector<RecurrenceTask> tasks;
};

#endif
