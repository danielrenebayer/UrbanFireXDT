#!/bin/bash
#
# measures the time required for compiling
# the complete simulation with optimization
#

main () {
    make opti
}

cd $(dirname $0)
cd ../code

time main

