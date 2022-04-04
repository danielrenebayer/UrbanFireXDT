#!/bin/bash
#
# measures the time required for compiling
# the complete simulation
#

main () {
    make
}

cd $(dirname $0)
cd ../code

time main

