#!/bin/bash

if [[ $EUID -ne 0 ]]; then
   echo "This script has to be run as root" 
   echo "to clear the linux page cache with /proc/sys/vm/drop_caches"
   exit 1
fi

echo "Compiling main.cpp..."
g++ -O3 -std=c++11 main.cpp -o runexperiment

echo "Running experiment... (this could take a while)"
./runexperiment | tee experiment.log

echo "Creating CSV with the results (results.csv)..."
cat experiment.log | grep @ | sed -n "s/@//p" > results.csv

echo "All done!"
