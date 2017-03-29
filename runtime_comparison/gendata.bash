#!/bin/bash
echo "Compiling gendata.cpp..."
g++ -O3 -std=c++11 gendata.cpp -o gendata
echo "Generating database.txt..."
echo "This file takes 6.4GB on disk"
./gendata
echo "Done!"
