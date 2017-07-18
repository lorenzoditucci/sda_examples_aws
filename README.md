## sda_examples_aws
Examples for SDAccel 2017.1+ on AWS F1 instances for Coursera

# Smith-Waterman
Porting Smith-Waterman code to SDAccel 2017.1 on KU3 with 2 Memory Ports.

Current Performance: 42.3 GCUPS

Paper [Here](http://ieeexplore.ieee.org/abstract/document/7927082/)


## dual-channel-kintex-smithwaterman_sdx_prj
Smith-Waterman deployed with GUI

Curent Problem: The gui does not offer a way to map ports

# Vector Addition
## makefile_vector_add
Simple vector addition deployed with the Makefile Flow

## vector_add
Same vector addition as above, deployed with GUI

# ProFAX

Hardware Acceleration of a Protein Folding Algorithm. Winning Entry of the Xilinx Open Hardware 2016, more info [here] (http://www.openhw.eu/2016-finalists.html)

Official Repository [here](https://bitbucket.org/necst/profax-src)

Paper [here](http://ieeexplore.ieee.org/abstract/document/7740584/)

# Makefile
Custom Makefile, change the source files and options for Xocc to run it.

```C
make [emulation | build | clean | clean_sw_emu | clean_hw_emu | clean_hw | cleanall] TARGET=<sw_emu | hw_emu | hw>
 
```

# Contacts & Contributions

Contributions to this repo are more than welcome!

For any request shoot an email at lorenzo.ditucci@polimi.it
