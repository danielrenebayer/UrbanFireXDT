#!/bin/bash
#
# Use this script to run all test simulations
# and compare the results with the verified output
#

# change to the directory where this script is located
cd $(dirname $0)

dirs_to_compare="S0001 S0002 S0003"
# 0. delete output dirs if they already exist (e.g. from previous test runs)
for di in $dirs_to_compare; do
    rm -rf test-output/$di
done

# 1. run simulations
../code/main-opti -m --config ../test/test-config/test_config.json 1
../code/main-opti -m --config ../test/test-config/test_config.json 2
../code/main-opti -m --config ../test/test-config/test_config.json 3

echo -e "\n------------------------------"
echo -e "-- Simulation runs finished --"
echo -e "------------------------------\n"

# 2. compare results with verified output
error_happend=0
for di in $dirs_to_compare; do
    #
    diff -rq test-output/$di test-output-verified/$di
    #
    # check if dirs are the same, if yes, return value is 0
    if (( $? != 0 )); then
        echo "ERROR for scenario $di: outputs do not match!"
        error_happend=1
    fi
done

# 3. output results
echo -e "\n------------------------------------"
echo "-- Test result:                   --"
if (( $error_happend > 0 )); then
    echo -e "--\033[01;31m Failed! Outputs do not match!\033[0m  --"
else
    echo -e "--\033[01;32m Everything passed!           \033[0m  --"
fi
echo -e "------------------------------------\n"

