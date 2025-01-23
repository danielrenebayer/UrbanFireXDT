/*
 * worker_threads.hpp
 *
 * It contains the definition of classes required for
 * the worker threads.
 *
 */

#ifndef WORKER_THREADS_HPP
#define WORKER_THREADS_HPP

#include <atomic>
#include <condition_variable>
#include <latch>
#include <list>
#include <mutex>
#include <thread>
#include <vector>

// The following classes are defined in this header file:
class CUControllerThreadGroupManager;
class CUControllerWorkerThread;

#include "global.h"
#include "units.h"

/*!
 * This class represents one group of threads of the class CUControllerWorkerThread.
 * It provides functionality for creating and managing all existing thread instances.
 * 
 * This class is not thread-safe.
 */
class CUControllerThreadGroupManager {
    private:
        friend class CUControllerWorkerThread;

    public:
        /*!
         * This function initializes the a new instance of the thread group manager.
         * It initializes the threads internally.
         * 
         * @param subsection: If the thread group should only be started for a subset of the available control units, the list of control units must be passed by this argument. Defaults to NULL, meaning that the group manager should handle all control units.
         */
        CUControllerThreadGroupManager(const std::vector<ControlUnit*>* subsection = NULL);

        ~CUControllerThreadGroupManager();
        //static void InitializeThreadGroupManager();
        
        /*!
         * This function starts all worker threads.
         * This means that the threads are forked and start to idle and with for tasks.
         * If called multiple times, it does not do anything.
         */
        void startAllWorkerThreads();

        /*!
         * This function stops all worker threads.
         * It stops all threads and joins them again.
         * If called multiple times, it does not do anything.
         */
        void stopAllWorkerThreads();

        /*!
         * This method notifies the worker threads to start working.
         * This means that the threads will start to call ControlUnit::compute_next_value() iteratively for all control units connected.
         * 
         * @param ts: The current time step (required for calling ControlUnit::compute_next_value())
         * @param dayOfWeek_l: The day of the week (left aligned) of the current time step (required for calling ControlUnit::compute_next_value())
         * @param hourOfDay_l: The hour of the day (left aligned) of the current time step (required for calling ControlUnit::compute_next_value())
         */
        void executeOneStep(unsigned long ts/*, const std::vector<ControlUnit*>* subsection = NULL*/);
        
        /*!
         * This function waits until all workers are finished with executed their task,
         * that has been started with ExecuteOneStep().
         */
        void waitForWorkersToFinish();

        /*!
         * This function deletes all created worker threads
         */
        // static void Vacuum();

    private:
        //bool initialized; /* = false */ ///< True, if the manager is already initialized
        std::vector<CUControllerWorkerThread*> worker_threads; ///< Vector of worker threads
        //static std::barrier<>* all_workers_finished_barrier;
        std::latch* all_workers_finished_latch;
        //static std::condition_variable cv_finished_signaling; ///< Internal conditional variable, required for notifying that a worker thread is finished
};

/*!
 * This class represents a working thread that calls the method
 * ControlUnit::compute_next_value() for all connected units on request.
 * The thread sleeps until it is activated using the public method
 * executeOneStepForAllConnCUs().
 */
class CUControllerWorkerThread {
    public:
        /*!
         * Constructs a new working thread for a list of controlled control units (CUs).
         * Attention: A control unit MUST ONLY be connected to ONE working thread! Otherwise ControlUnit::compute_next_value()
         * is called more than once for the control units connected!
         * 
         * @param connected_units_: The list of control units that are connected to this working thread
         * @param caller: The reference to the thread group manager object
         */
        CUControllerWorkerThread(const std::list<ControlUnit*>* connected_units_, CUControllerThreadGroupManager& caller);
        ~CUControllerWorkerThread();
        void start(); ///< Starts this thread. This method has to be called before the call of executeOneStepForAllConnCUs().
        void stop(); ///< Stops this working thread.

        /*!
         * This method notifies the thread to start working.
         * This means that the thread will start to call ControlUnit::compute_next_value() iteratively for all control units connected.
         * Basically, it sets atomic_flag_exec to true.
         * 
         * @param ts: The current time step (required for calling ControlUnit::compute_next_value())
         * @param dayOfWeek_l: The day of the week (left aligned) of the current time step (required for calling ControlUnit::compute_next_value())
         * @param hourOfDay_l: The hour of the day (left aligned) of the current time step (required for calling ControlUnit::compute_next_value())
         * @param subsection:  If the simulation should only be executed for a subset of instances, it can be defined here
         */
        void executeOneStepForAllConnCUs(unsigned long ts/*, const std::vector<ControlUnit*>* subsection = NULL*/);
        
        /*!
         * This function returns true if the thread is idling, i.e., it is running, but not currently working and has also no planned work.
         */
        bool isIdling() const { return atomic_flag_idling; }

    private:
        CUControllerThreadGroupManager& thread_group_manager; ///< Internal reference to the thread group manager
        std::vector<ControlUnit*> connected_units; ///< List of connected control units
        std::thread current_thread;
        std::mutex mtx; ///< The mutex object per instance
        std::condition_variable cv; ///< The conditional variable to signal the requests for running, i.e., by calling executeOneStepForAllConnCUs()
        std::atomic<bool> atomic_flag_stop;    ///< Set to true if the the working thread should stop
        std::atomic<bool> atomic_flag_exec;    ///< Set to true if the main action of this worker (i.e., calling ControlUnit::compute_next_value() ) should be executed
        std::atomic<bool> atomic_flag_running; ///< Set to true if the thread is invoked and running (working or idling)
        std::atomic<bool> atomic_flag_idling;  ///< Set to true if the thread is idling (i.e., running but without work and without planned work)
        //
        std::atomic<unsigned long> atomic_param_ts; ///< Parameter required for calling ControlUnit::compute_next_value()
        //
        void run(); ///< Main internal function for this thread. It is started and stopped with the start()- and stop()-method.
};

#endif

