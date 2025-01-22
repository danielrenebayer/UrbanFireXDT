/**
 * optimization_unit_general.hpp
 *
 * This file contains all general classes / structs required
 * by the optimization implementations.
 */

#ifndef OPTIMIZATION_UNIT_GENERAL_HPP
#define OPTIMIZATION_UNIT_GENERAL_HPP

#include <vector>

namespace CUOptimization {

    class OptimalControlCommandsOverHorizon {
        public:
            std::vector<double> hp_power; ///< The heat pump power for every time step in the time horizon
            std::vector<double> cs_power; ///< The charging station power for every time step in the time horizon
            std::vector<double> bs_power; ///< The battery storage power for every time step in the time horizon

            /**
             * Constructs a new object with 0.0-initialized vectors with length of parameter time_horizon
             */
            OptimalControlCommandsOverHorizon(unsigned int time_horizon) :
                hp_power(time_horizon, 0.0), cs_power(time_horizon, 0.0), bs_power(time_horizon, 0.0)
                {}

            /**
             * Resets the vectors with 0.0-initialized values with length of parameter (new) parameter time_horizon
             */
            void reset(unsigned int new_horizon) {
                hp_power.assign(new_horizon, 0.0);
                cs_power.assign(new_horizon, 0.0);
                bs_power.assign(new_horizon, 0.0);
            }

            /**
             * Removes the first element, shifts the second one to the first place and so on and
             * adds a 0.0 at the last place.
             */
            void shiftVectorsByOnePlace() {
                std::move(hp_power.begin() + 1, hp_power.end(), hp_power.begin());
                hp_power.back() = 0.0;
                std::move(cs_power.begin() + 1, cs_power.end(), cs_power.begin());
                cs_power.back() = 0.0;
                std::move(bs_power.begin() + 1, bs_power.end(), bs_power.begin());
                bs_power.back() = 0.0;
            }
    };
}

#endif

