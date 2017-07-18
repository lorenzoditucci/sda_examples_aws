#BASH shell setup script for compiling application using xocc and running
#CPU and hardware emulation directly from command line
#
#set XILINX_OPENCL to SDACCel installation directory
export XILINX_OPENCL=/opt/Xilinx/SDAccel/2015.4

source $XILINX_OPENCL/settings64.sh
export LD_LIBRARY_PATH=$XILINX_OPENCL/runtime/lib/x86_64:$LD_LIBRARY_PATH
export PATH=$XILINX_OPENCL/lnx64/tools/gcc/bin/:$PATH
