#!/bin/bash

#
# Use this file to create the doxygen documentation locally
#

# change directory to this folder, wherever it is called from
cd $(dirname $0)

doxygen doxygen.conf

