#!/bin/bash
#
# Use this script to run all test simulations
# and compare the results with the verified output
#

# change to the directory where this script is located
cd $(dirname $0)

dirs_to_compare="S0001 S0002 S0003 S0004 S0005 S0006 S0007 S0008"
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
../bin/simulation-opti --config ../test/test-config/test_config.json 1
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 2. Simulation run --"
../bin/simulation-opti --config ../test/test-config/test_config.json 2
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 3. Simulation run --"
../bin/simulation-opti --config ../test/test-config/test_config.json 3
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 4. Simulation run --"
../bin/simulation-opti --config ../test/test-config/test_config.json 4
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 5. Simulation run --"
../bin/simulation-opti --config ../test/test-config/test_config.json 5
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 6. Simulation run --"
../bin/simulation-opti --config ../test/test-config/test_config.json --pvar 1 7
if (( $? != 0 )); then program_error=1; fi
echo -e "\n-- 7. Simulation run --"
../bin/simulation-opti --config ../test/test-config/test_config.json --seed 1234 --ev-output all 8
if (( $? != 0 )); then program_error=1; fi


# use valgrind for last run, but only, if there are no errors before
if (( $program_error <= 0 )); then
    echo -e "\n-- 8. Simulation run (using valgrind) --"
    valgrind --error-exitcode=1 ../bin/simulation-dbg --config ../test/test-config/test_config.json --seed 1234 6
    if (( $? != 0 )); then memory_error=1; else memory_error=0;  fi
fi

echo -e "\n------------------------------"
echo -e "-- Simulation runs finished --"
echo -e "------------------------------\n"

# 2. compare results with verified output
output_error=0
for di in $dirs_to_compare; do
    #
    # Sort the file ST1-AllCUs-ts.csv, as the output can be mixed up due to parallelization
    find test-output/$di -name ST1-AllCUs-ts.csv -type f -exec sh -c '
for file; do
    sort "$file" > "${file}-sorted"
    done
' sh {} +
    # Sort ev-details.csv
    find test-output/$di -name ev-details.csv -type f -exec sh -c '
for file; do
    sort "$file" > "${file}-sorted"
    done
' sh {} +
    #
    # Do the main diff task (but ignore ST1-AllCUs-ts.csv)
    diff -rq --exclude=build_and_run_info.txt --exclude=runtime-information.csv --exclude=ST1-AllCUs-ts.csv --exclude=ev-details.csv test-output/$di test-output-verified/$di
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
