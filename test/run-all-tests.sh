#!/bin/bash
#
# Use this script to run all test simulations
# and compare the results with the verified output
#

# change to the directory where this script is located
cd $(dirname $0)

dirs_to_compare="S0001 S0002 S0003 S0004"
# 0. delete output dirs if they already exist (e.g. from previous test runs)
for di in $dirs_to_compare; do
    rm -rf test-output/$di
done

# 0b. Define the errors that can happen
program_error=0  # An error during the simulation execution
memory_error=-1  # An error detected by valgrind
output_error=0   # Outputs do not match the verified output

# 1. run simulations
echo -e "\n-- 1. Simulation run --"
../code/main-opti --config ../test/test-config/test_config.json 1
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 2. Simulation run --"
../code/main-opti --config ../test/test-config/test_config.json 2
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 3. Simulation run --"
../code/main-opti --config ../test/test-config/test_config.json 3
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 4. Simulation run --"
../code/main-opti --config ../test/test-config/test_config.json 4
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 5. Simulation run --"
../code/main-opti --config ../test/test-config/test_config.json 5
if (( $? != 0 )); then program_error=1; fi


# use valgrind for last run, but only, if there are no errors before
if (( $program_error <= 0 )); then
    echo -e "\n-- 6. Simulation run (using valgrind) --"
    valgrind ../code/main-dbg --config ../test/test-config/test_config.json 6
    if (( $? != 0 )); then memory_error=1; else memory_error=0;  fi
fi

echo -e "\n------------------------------"
echo -e "-- Simulation runs finished --"
echo -e "------------------------------\n"

# 2. compare results with verified output
output_error=0
for di in $dirs_to_compare; do
    #
    diff -rq --exclude=build_and_run_info.txt test-output/$di test-output-verified/$di
    #
    # check if dirs are the same, if yes, return value is 0
    if (( $? != 0 )); then
        echo -e "\n\033[01;31mERROR for scenario $di: outputs do not match!\033[0m\n"
        output_error=1
    fi
done

# 3. output results
echo -e "\n--------------------------------------"
echo "-- Simulation procedure:            --"
if (( $program_error > 0 )); then
    echo -e "--\033[01;31m One or more errors occured     \033[0m  --"
    echo -e "--\033[01;31m during the simulation run.     \033[0m  --"
else
    echo -e "--\033[01;32m Works without errors!          \033[0m  --"
fi
echo -e "--------------------------------------"
echo "-- Output test result:              --"
if (( $output_error > 0 )); then
    echo -e "--\033[01;31m Failed! Outputs do not match!  \033[0m  --"
else
    echo -e "--\033[01;32m Everything passed!             \033[0m  --"
fi
echo -e "--------------------------------------"
echo "-- Memory check resuly:             --"
if (( $program_error > 0 )); then
  echo -e "-- not evaluated, as errors occured above   --"
else
  if   (( $memory_error < 0 )); then
    echo -e "--\033[01;33m Not checked due to prev. errors! \033[0m--"
  elif (( $memory_error > 0 )); then
    echo -e "--\033[01;31m Failed! Memory lekage detected!  \033[0m--"
  else
    echo -e "--\033[01;32m Everything passed!             \033[0m  --"
  fi
fi
echo -e "--------------------------------------\n"
