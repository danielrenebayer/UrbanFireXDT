import os
import sys
import time

sys.path.insert(0, os.path.abspath("../bin"))

import UrbanFireXDT

# Initialize the C++ simulation
UrbanFireXDT.initialize({"config":"test-config/test_config.json", "scenario": 1})

# The number of episodes (i.e., repetitions of the complete simulation time span) to execute
n_episodes = 3
simulation_needs_reset = False

# Loop over all episodes
for episode_id in range(1, n_episodes+1):
    print(f"Starting new episode {episode_id} ...")

    # Reset simulation to the starting state, if a new episode is starting (not required for the first time)
    if simulation_needs_reset:
        UrbanFireXDT.reset_simulation_to_start()
    simulation_needs_reset = True # Reset in future runs

    # Loop over all steps in the simulation until it reaches the end
    while not UrbanFireXDT.simulation_finished():
        this_timestep_id = UrbanFireXDT.get_next_timestep_id()
        print(f"Simulation step {this_timestep_id} for episode {episode_id} ...", end="")

        # 1. Query the state from C++
        cuID = 1
        current_state = UrbanFireXDT.get_state(controlUnitID=cuID)
        #print(f"C++ state of control unit with ID {cuID} = {current_state}")

        # 2. Determine and send new commands
        new_commands = {"test_cmd": 1234.0}
        UrbanFireXDT.send_commands(controlUnitID=cuID, commands=new_commands)

        # 3. Tell C++ to run the next step
        UrbanFireXDT.run_one_step()
        print(" -> .")

    print(f"Episode {episode_id} is finished.")

print("Simulation clean up started ...")
UrbanFireXDT.vacuum()
print("Simulation finished and vacuumed.")

