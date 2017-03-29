Tool to compare quality of SSJ, WS-join, US-join, HWS-join and HSSJ.

To compile and run the experiments, you can use `runexperiments.bash`:
```bash
./runexperiments.bash
```

Make sure that enough memory is available on your machine! Approximately 11 * n1 * 64 bits of memory are needed to run the experiments, for the default value of n1 this corresponds to 17.6 GB of memory. If desired, experiment parameters can be changed directly in `qualityComparison.cpp`.

`qualityComparison.cpp`
> compare estimation quality 

`sampleJoins.h`
> join sampling algorithms and utilities

`mtwist.h`
> header-only implementation of a mersenne prime twister

`./runexperiments.bash`
> compile and run quality experiments
