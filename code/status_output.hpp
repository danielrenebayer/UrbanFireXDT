/*
 * status_output.hpp
 *
 * This file contains all code for the status output of the simulation
 * using the ncurses library.
 *
 */

#ifndef STATUS_OUTPUT_HPP
#define STATUS_OUTPUT_HPP

#include <ncurses.h>
#include <string>
#include <vector>
#include <mutex>
#include <thread>
#include <atomic>
#include <fstream>

/**
 * @class StatusOutput
 * @brief A thread-safe logging and status output system using ncurses.
 *
 * This class provides a singleton-style interface for handling real-time
 * status output in a terminal window divided into three panes:
 * - Left: General simulation output.
 * - Top-right: Error messages.
 * - Bottom-right: Current status updated periodically.
 *
 * It also logs messages to separate log files.
 */
class StatusOutput {
    public:
        /**
         * @brief Adds a message to the simulation output log.
         * @param message The message to be logged.
         */
        static void add_status_output(const std::string& message);

        /**
         * @brief Adds an error message to the error log.
         * @param error The error message to be logged.
         */
        static void add_error_message(const std::string& error);

        /**
         * @brief Initializes the ncurses interface.
         *
         * This must be called at the start of main() before any output is displayed.
         */
        static void initialize_ncurses();

        /**
         * @brief Shuts down the ncurses interface.
         *
         * This should be called before the program terminates to restore the terminal state.
         */
        static void shutdown_ncurses();

        /**
         * @brief Starts the background thread that periodically updates the status display.
         *
         * This should be called when the main simulation run starts.
         */
        static void start_status_updater_thread();

        /**
         * @brief Stops the background status updater thread.
         *
         * This should be called when the simulation ends to properly clean up resources.
         */
        static void stop_status_updater_thread();

    private:
        static void update_status_window();
        static void status_updater();

        static std::mutex log_mutex;
        static std::mutex error_mutex;
        static std::mutex status_mutex;

        static std::vector<std::string> log_messages;
        static std::vector<std::string> error_messages;
        static std::string current_status;

        static std::atomic<bool> running;
        static std::atomic<bool> ncurses_initialized;
        static std::thread updater_thread;

        static WINDOW* log_win;
        static WINDOW* error_win;
        static WINDOW* status_win;
};

#endif // STATUS_OUTPUT_HPP
