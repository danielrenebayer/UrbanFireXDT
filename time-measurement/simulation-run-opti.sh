#!/bin/bash
#
# measures the time required for one complete
# run of the compiler-optimized simulation
#

main () {
    ./main-opti
}

cd $(dirname $0)
cd ../code

time main

