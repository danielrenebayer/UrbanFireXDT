#!/bin/bash
#
# measures the time required for compiling
# the complete simulation with debug marks
#

main () {
    make debug
}

cd $(dirname $0)
cd ../code

time main

