Tool to compare runtime of SSJ, WS-join, US-join, HWS-join and HSSJ.
This tool is Linux only as it depends on system calls to `/proc/sys/vm/drop_caches` (among others).

To run experiments, you have to first generate data and then run the experiments. You need root to run the experiments, since these rights are required to drop the page cache.

```bash
./gendata.bash
sudo ./runexperiments.bash
```
Note that ample RAM is needed to store all columns in memory and that the flushing code assumes that the L3 cache is much smaller than 50MiB. If desired, experiment parameters can be changed directly in `main.cpp`. 

`main.cpp`
> code to run benchmarks

`gendata.cpp`
> can be used to generate data on disk

`runexperiments.bash`
> script that compiles and runs benchmarks and creates a CSV file with the results

`gendata.bash`
> script that compiles and runs `gendata.cpp` to generate a database (requires 6.4 GB of diskspace)


