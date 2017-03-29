#!/bin/bash
echo "compiling qualityComparison.cpp ..."
g++ -O3 -std=c++11 qualityComparison.cpp -o qualityComparison
echo "running experiments ..."
./qualityComparison
