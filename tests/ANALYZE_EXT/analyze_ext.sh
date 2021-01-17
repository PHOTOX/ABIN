#!/bin/bash

# This is an external script that is called every time step from analyze_ext.F90.
# You can use it e.g. for on-the-fly analysis of ab initio output.

# Here we just test that it is being called.
# Currently, we pass the time step as a BASH parameter from ABIN
# so that we can decide how often to perform the analysis.
timestep=$1

echo "Hello from analyze_ext.sh! Time step=$timestep" >> analyze_ext.dat

exit 0
