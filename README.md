## sda_examples_aws
Examples for SDAccel 2017.1+ on AWS F1 instances

ap_int problem fixed by changing compiler in the Makefile as in the SDAccel Eamples repo

# Makefile
Custom Makefile, change the source files and options for Xocc to run it.

```C
make [emulation | build | clean | clean_sw_emu | clean_hw_emu | clean_hw | cleanall] TARGET=<sw_emu | hw_emu | hw>
 
```

