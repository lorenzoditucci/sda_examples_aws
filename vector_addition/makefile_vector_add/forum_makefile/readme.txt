This is an exmaple showing how to use Makefile and xocc to compile OpenCL host
code using gcc and kernel code using xocc.

The example also includes environment setup scripts setup_xocc.csh (Csh/Tcsh)
and setup_xocc.sh (Bash) for running CPU and hardware emulation. Note that 
the included setup scripts do not work for running application on hardware. 
To run application on hardware, the setup script generated from Xilinx board
installation utility "xbinst" must be sourced.

CPU and Hardware emulation
-------------------------------------------------
* Set up environment
$source setup_xocc.csh
or
$source setup_xocc.sh

* Compile and run CPU emulation
$make

* Compile and run hardware emulation
$make FLOW=hw_emu

Running on Hardware
-------------------------------------------------
* Set up environment
$source <xbinst_dir>/setup_xocc.csh
or
$source <xbinst_dir>/setup_xocc.sh

* Compile and run on hardware
$make FLOW=hw


