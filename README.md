## sda_examples_aws
Examples for SDAccel 2017.1+ on AWS F1 instances

# Smith-Waterman
Porting Smith-Waterman code to SDAccel 2017.1 on KU3 with 2 Memory Ports.

Current Performance: 42.3 GCUPS

Paper [Here](http://ieeexplore.ieee.org/abstract/document/7927082/)

@INPROCEEDINGS{7927082, 
author={L. Di Tucci and K. O'Brien and M. Blott and M. D. Santambrogio}, 
booktitle={Design, Automation Test in Europe Conference Exhibition (DATE), 2017}, 
title={Architectural optimizations for high performance and energy efficient Smith-Waterman implementation on FPGAs using OpenCL}, 
year={2017}, 
pages={716-721}, 
keywords={Algorithm design and analysis;Computer architecture;Field programmable gate arrays;Graphics processing units;Hardware;Optimization;Performance evaluation}, 
doi={10.23919/DATE.2017.7927082}, 
month={March},}

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

Citing us:

@INPROCEEDINGS{7740584, 
author={G. Guidi and L. Di Tucci and M. D. Santambrogio}, 
booktitle={2016 IEEE 2nd International Forum on Research and Technologies for Society and Industry Leveraging a better tomorrow (RTSI)}, 
title={ProFAX: A hardware acceleration of a protein folding algorithm}, 
year={2016}, 
pages={1-6}, 
keywords={Monte Carlo methods;drugs;macromolecules;medical computing;molecular biophysics;proteins;Monte Carlo simulation;ProFAX;ab initio protein folding algorithm;amino acid sequence;drug industries;energetic features;geometrical features;hardware acceleration;pharmaceutical therapies;software implementation;tertiary structure;Acceleration;Algorithm design and analysis;Amino acids;Field programmable gate arrays;Hardware;Industries;Proteins}, 
doi={10.1109/RTSI.2016.7740584}, 
month={Sept},}

# Makefile
Custom Makefile, change the source files and options for Xocc to run it.

```C
make [emulation | build | clean | clean_sw_emu | clean_hw_emu | clean_hw | cleanall] TARGET=<sw_emu | hw_emu | hw>
 
```

# Contacts & Contributions

Contributions to this repo are more than welcome!

For any request shoot an email at lorenzo.ditucci@polimi.it
