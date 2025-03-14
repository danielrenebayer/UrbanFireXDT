#include "status_output.hpp"

#include <iostream>
#include <chrono>
#include <sstream>

#include "global.h"
#include "output.h"

// Define static members
std::mutex StatusOutput::log_mutex;
std::mutex StatusOutput::error_mutex;
std::mutex StatusOutput::status_mutex;

std::vector<std::string> StatusOutput::log_messages;
std::vector<std::string> StatusOutput::error_messages;
std::string StatusOutput::current_status;

std::atomic<bool> StatusOutput::running{false};
std::atomic<bool> StatusOutput::ncurses_initialized{false};
std::thread StatusOutput::updater_thread;

WINDOW* StatusOutput::log_win = nullptr;
WINDOW* StatusOutput::error_win = nullptr;
WINDOW* StatusOutput::status_win = nullptr;

// Main implementation part
void StatusOutput::initialize_ncurses() {
    if (!ncurses_initialized) {
        initscr();
        noecho();
        cbreak();
        refresh();
        
        int height, width;
        getmaxyx(stdscr, height, width);
        
        log_win = newwin(height, width / 2, 0, 0);
        error_win = newwin(height / 2, width / 2, 0, width / 2);
        status_win = newwin(height / 2, width / 2, height / 2, width / 2);
        
        box(log_win, 0, 0);
        box(error_win, 0, 0);
        box(status_win, 0, 0);
        
        wrefresh(log_win);
        wrefresh(error_win);
        wrefresh(status_win);

        ncurses_initialized = true;
    }
}

void StatusOutput::shutdown_ncurses() {
    if (ncurses_initialized) {
        endwin();
        ncurses_initialized = false;
    }
}

void StatusOutput::add_status_output(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    log_messages.push_back(message);
    /*
    if (output::logstream_status != NULL) {
       *output::logstream_status << message << "\n";
        output::logstream_status->flush();
    }
    */
}

void StatusOutput::add_error_message(const std::string& error) {
    std::lock_guard<std::mutex> lock(error_mutex);
    error_messages.push_back(error);
    /*
    if (output::logstream_errors != NULL) {
       *output::logstream_errors << error << "\n";
        output::logstream_errors->flush();
    }
    */
}

void StatusOutput::update_status_window() {
    std::lock_guard<std::mutex> lock(status_mutex);
    werase(status_win);
    box(status_win, 0, 0);
    mvwprintw(status_win, 1, 1, "Status: %s", current_status.c_str());
    wrefresh(status_win);
}

void StatusOutput::status_updater() {
    while (running) {
        {
            auto time_now = std::chrono::system_clock::now();
            auto time_diff = std::chrono::duration_cast<std::chrono::seconds>(time_now - global::time_of_simulation_start).count();
            //
            std::lock_guard<std::mutex> lock(status_mutex);
            current_status  = "Time since simulation start = ";
            current_status += to_string( time_diff );
            current_status += "s - ";
            current_status += "Number of optimization calls = ";
            current_status += to_string( ControlUnit::GetTotalNumberOfOptimizationCalls() );
            if (Global::get_max_parallel_opti_vars() > 0) {
                current_status += " with ";
                current_status += to_string( ControlUnit::GetCurrentNumberOfOptiVars() );
                current_status += " parallel optimization variables";
            }
            current_status += "\r";
            //
            ofstream ofs("/tmp/simulation-status.txt", std::ofstream::out);
            ofs << "Time since simulation start  = " << to_string( time_diff ) << "s\n";
            ofs << "Number of optimization calls = " << to_string( ControlUnit::GetTotalNumberOfOptimizationCalls() ) << "\n";
            if (Global::get_max_parallel_opti_vars() > 0) {
                ofs << "Number of parallel optimization vars     = " << to_string( ControlUnit::GetCurrentNumberOfOptiVars() ) << "\n";
            }
            ofs << "ControlUnit::GetNumberOfCUsWithSimCompPV = " << to_string( ControlUnit::GetNumberOfCUsWithSimCompPV()  ) << "\n";
            ofs << "ControlUnit::GetNumberOfCUsWithSimCompHP = " << to_string (ControlUnit::GetNumberOfCUsWithSimCompHP() ) << "\n";
            ofs << "ControlUnit::GetNumberOfCUsWithSimCompEV = " << to_string( ControlUnit::GetNumberOfCUsWithSimCompEV() ) << "\n";
            ofs.close();
        }
        if (ncurses_initialized) {
            update_status_window();
        } else {
            // writing to stdout
            std::cout << current_status << std::flush;
        }
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }
}

void StatusOutput::start_status_updater_thread() {
    if (!running) {
        running = true;
        updater_thread = std::thread(status_updater);
    }
}

void StatusOutput::stop_status_updater_thread() {
    if (running) {
        running = false;
        updater_thread.join();
    }
}
