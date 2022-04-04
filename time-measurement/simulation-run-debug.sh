#!/bin/bash
#
# measures the time required for one complete
# run of the simulation with debug marks
#

main () {
    ./main-dbg
}

cd $(dirname $0)
cd ../code

time main

