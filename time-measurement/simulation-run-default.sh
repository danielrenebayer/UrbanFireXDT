#!/bin/bash
#
# measures the time required for one complete
# run of the simulation
#

main () {
    ./main
}

cd $(dirname $0)
cd ../code

time main

