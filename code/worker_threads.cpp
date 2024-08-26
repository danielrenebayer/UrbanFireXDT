#include "worker_threads.hpp"

#include <iostream>
#include <latch>
#include <list>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include "global.h"
#include "units.h"


// ----------------------------------- //
//         Implementation of           //
//   CUControllerThreadGroupManager    //
// ----------------------------------- //
bool CUControllerThreadGroupManager::initialized = false;
std::vector<CUControllerWorkerThread*> CUControllerThreadGroupManager::worker_threads;
std::latch* CUControllerThreadGroupManager::all_workers_finished_latch = NULL;

void CUControllerThreadGroupManager::InitializeThreadGroupManager() {
    if (!initialized) {
        initialized = true;

        if (Global::get_n_threads() < 1)
            return;

        //all_workers_finished_barrier = new std::barrier<>(Global::get_n_threads() + 1); // +1 as the main thread calls arrive_and_wait(), too

        // create as much lists as we have workers where we iteratively add one control unit
        std::vector<std::list<ControlUnit*>> vlc;
        for (unsigned int t = 0; t < Global::get_n_threads(); t++) {
            vlc.emplace_back(); // create an empty list per worker
        }
        // add control units to the lists iteratively
        unsigned int t = 0; // the current worker where the next control unti will be attached to
        for (ControlUnit* cu : ControlUnit::GetArrayOfInstances()) {
            vlc[t].push_back(cu);
            // increment t
            t++;
            if ( t >= Global::get_n_threads() )
                t = 0;
        }

        worker_threads.reserve( Global::get_n_threads() );
        for (unsigned int t = 0; t < Global::get_n_threads(); t++) {
            CUControllerWorkerThread* new_thread = new CUControllerWorkerThread( &vlc[t] );
            worker_threads.push_back(new_thread);
        }
    }
}

void CUControllerThreadGroupManager::StartAllWorkerThreads() {
    for (CUControllerWorkerThread* wt : worker_threads) {
        wt->start();
    }
}

void CUControllerThreadGroupManager::StopAllWorkerThreads() {
    for (CUControllerWorkerThread* wt : worker_threads) {
        wt->stop();
    }
}

void CUControllerThreadGroupManager::ExecuteOneStep(unsigned long ts, unsigned int dayOfWeek_l, unsigned int hourOfDay_l) {
    // (Re-)Initialize the latch object
    if (all_workers_finished_latch != NULL)
        delete all_workers_finished_latch;
    all_workers_finished_latch = new std::latch(Global::get_n_threads());
    // Notify all threads to start working
    for (CUControllerWorkerThread* wt : worker_threads) {
        wt->executeOneStepForAllConnCUs(ts, dayOfWeek_l, hourOfDay_l);
    }
}

void CUControllerThreadGroupManager::WaitForWorkersToFinish() {
    // for (CUControllerWorkerThread* wt : worker_threads)
    //all_workers_finished_barrier->arrive_and_wait();
    all_workers_finished_latch->wait();
}

void CUControllerThreadGroupManager::Vacuum() {
    if (initialized) {
        for (CUControllerWorkerThread* thread : worker_threads) {
            delete thread;
        }
        worker_threads.clear();
        initialized = false;
        //delete all_workers_finished_barrier;
        delete all_workers_finished_latch;
    }
}



// ----------------------------- //
//      Implementation of        //
//   CUControllerWorkerThread    //
// ----------------------------- //
CUControllerWorkerThread::CUControllerWorkerThread(std::list<ControlUnit*>* connected_units_) {
    // add the stations to connect to this thread to the internal vector storing all stations
    // check, if a station is not already connected to some other thread
    connected_units = std::vector<ControlUnit*>();
    for (ControlUnit* s : *connected_units_) {
        if (s->worker_thread == NULL) {
            s->worker_thread = this;
            connected_units.push_back(s);
        } else {
            std::cerr << "Error: Control Unit with ID " << std::to_string(s->get_unitID()) << " already connected to another working thread!" << std::endl;
        }
    }
    //
    atomic_flag_stop    = false;
    atomic_flag_exec    = false;
    atomic_flag_running = false;
    atomic_flag_idling  = true;
}

CUControllerWorkerThread::~CUControllerWorkerThread() {
    // Stop the thread and clean up
    stop();
    if (current_thread.joinable()) {
        current_thread.join();
    }
    {
        std::unique_lock<std::mutex> lock_obj(mtx);
        atomic_flag_running = false;
    }
}

void CUControllerWorkerThread::start() {
    // Is the thread already started?
    bool start_thread = false;
    {
        std::unique_lock<std::mutex> lock_obj(mtx);
        if (!atomic_flag_running) {
            start_thread = true;
            atomic_flag_running = true;
        }
    }
    // Start the thread, if it is not already running
    if (start_thread) {
        current_thread = std::thread(&CUControllerWorkerThread::run, this);
    }
}

void CUControllerWorkerThread::stop() {
    {
        // lock the mutex in this scope
        std::unique_lock<std::mutex> lock_obj(mtx);
        // do the atomic operation
        atomic_flag_stop = true;
    }
    cv.notify_all();
}

void CUControllerWorkerThread::executeOneStepForAllConnCUs(
    unsigned long ts,
    unsigned int dayOfWeek_l,
    unsigned int hourOfDay_l /*,
    const std::vector<ControlUnit*>* subsection *//* = NULL */)
{
    {
        // lock the mutex in this scope
        std::unique_lock<std::mutex> lock_obj(mtx);
        // set the flag
        atomic_flag_exec   = true;
        atomic_flag_idling = false;
        // set the variables storing the parameters for calling ControlUnit::compute_next_value()
        atomic_param_ts          = ts;
        atomic_param_dayOfWeek_l = dayOfWeek_l;
        atomic_param_hourOfDay_l = hourOfDay_l;
    }
    cv.notify_all();
}

void CUControllerWorkerThread::run() {
    bool exec_task = false; // thread-internal variable to store the value of atomic_flag_exec
    //
    unsigned long ts;          // local copies of the atomic variables required for passing them as an argument to ControlUnit::compute_next_value()
    unsigned int  dayOfWeek_l;
    unsigned int  hourOfDay_l;
    //
    while(true)
    {
        {
            // lock the mutex
            std::unique_lock<std::mutex> lock_obj(mtx);
            // wait for the variables to change ( the wait-method relases the mutex until the condition is met, otherwise other thrads could not aquire a lock in the executeOneStepForAllConnCUs()-method )
            cv.wait(lock_obj, [this] { return atomic_flag_stop || atomic_flag_exec; });

            // execute the main task if it is selected
            if (atomic_flag_exec) {
                exec_task = true;
                atomic_flag_exec = false;
                // copy the atomic variables to local copies
                ts          = atomic_param_ts;
                dayOfWeek_l = atomic_param_dayOfWeek_l;
                hourOfDay_l = atomic_param_hourOfDay_l;
            }

            // return if stop flag has been set
            if (atomic_flag_stop) {
                atomic_flag_stop = false;
                return;
            }
        }
        // call ControlUnit::compute_next_value() for all connected CUs
        // outside of the locked mutex
        if (exec_task) {
            exec_task = false; // do not execute it again
            //
            // iterate over all connected CUs
            // and execute the compute_next_value()-method
            for (ControlUnit* cu : connected_units) {
                cu->compute_next_value(ts, dayOfWeek_l, hourOfDay_l);
            }
        }
        // set atomic flag for idling to true AFTER the task has been executed
        {
            std::unique_lock<std::mutex> lock_obj(mtx);
            atomic_flag_idling = true;
        }
        //CUControllerThreadGroupManager::cv_finished_signaling.notify_all();
        // decrement the barrier
        //CUControllerThreadGroupManager::all_workers_finished_barrier->arrive();
        CUControllerThreadGroupManager::all_workers_finished_latch->count_down();
    }
}
