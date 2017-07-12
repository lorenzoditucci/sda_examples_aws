#######################################################################################################################################
#
#	Basic Makefile for SDAccel 2017.1
#	Usage make [emulation | build | clean | clean_sw_emu | clean_hw_emu | clean_hw | cleanall] target=<sw_emu | hw_emu | hw>
#
#
#######################################################################################################################################
XOCC=xocc
CC=g++

#Host code
host_src=./test-cl-2ddr.cpp
host_hdrs=
host_cflags=-D FPGA_DEVICE -g -Wall -I${XILINX_SDX}/runtime/include/1_2 -D C_KERNEL -O3 -Wall
host_lflags=-L${XILINX_SDX}/runtime/lib/x86_64 -lxilinxopencl

#name of host executable
host_exe=host_sw

#kernel
kernel_src=./smithwaterman.cpp
kernel_hdrs=
kernel_flags=
kernel_exe=smithwaterman
kernel_name=smithwaterman

#custom flag to give to xocc
kernel_LDCLFLAGS=--nk $(kernel_name):1 \
	--xp param:compiler.preserveHlsOutput=1 \
	--max_memory_ports $(kernel_name) \
	--memory_port_data_width $(kernel_name):512 \
	--xp misc:map_connect=add.kernel.smithwaterman_1.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI \
	--xp misc:map_connect=add.kernel.smithwaterman_1.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI \
	--xp misc:map_connect=add.kernel.smithwaterman_1.M_AXI_GMEM2.core.OCL_REGION_0.M00_AXI 


target_device=xilinx:aws-vu9p-f1:4ddr-xpr-2pr:4.0

#target for compilation [sw_emu | hw_emu | hw]
target=none
report_flag=n
report=
ifeq (${target}, sw_emu)
$(info software emulation)
target=sw_emu
ifeq (${report_flag}, y)
$(info creating report for software emulation set to true. This is going to take longer at it will synthesize the kernel)
report=--report estimate
else
$(info I am not creating a report for software emulation, set report_flag=y if you want it)
report=
endif
else ifeq (${target}, hw_emu)
$(info hardware emulation)
target=hw_emu
report=--report estimate
else ifeq (${target}, hw)
$(info system build)
target=hw
report=--report system
else
$(info no target selected)
endif


ifndef XILINX_SDX
$(error XILINX_SDX is not set. Please source the SDx settings64.{csh,sh} first)
endif

clean:
	rm -rf .Xil emconfig.json 

clean_sw_emu: clean
	rm -rf sw_emu
clean_hw_emu: clean
	rm -rf hw_emu
clean_hw: clean
	rm -rf hw

cleanall: clean_sw_emu clean_hw_emu clean_hw
	rm -rf _xocc_*

check_target:
ifeq (${target}, none)
	$(error Target can not be set to none)
endif

host:  check_target $(host_src) $(host_hdrs)
	mkdir -p $(target)
	$(CC) $(host_src) $(host_hdrs) $(host_cflags) $(host_lflags) -o $(target)/$(host_exe)

xo:	check_target
	mkdir -p $(target)
	$(XOCC) --platform $(target_device) --target $(target) --compile --include $(kernel_hdrs) --save-temps $(report) --kernel $(kernel_name) $(kernel_src) $(kernel_LDCLFLAGS) $(kernel_flags)  --output $(target)/$(kernel_exe).xo

xclbin:  check_target xo
	$(XOCC) --platform $(target_device) --target $(target) --link --include $(kernel_hdrs) --save-temps $(report) --kernel $(kernel_name) $(target)/$(kernel_exe).xo $(kernel_LDCLFLAGS) $(kernel_flags) --output $(target)/$(kernel_exe).xclbin

emulation:  host xclbin
	export XCL_EMULATION_MODE=$(target)
	emconfigutil --xdevice $(target_device) --nd 1
	./$(target)/$(host_exe) $(target)/$(kernel_exe).xclbin
	$(info Remeber to export XCL_EMULATION_MODE=$(target) and run emcondigutil for emulation purposes)

#sw_emu: host xclbin
#	export XCL_EMULATION_MODE=$(target)
#	emconfigutil --xdevice $(target_device) --nd 1
#	#emconfigutil -f xilinx:adm-pcie-ku3:2ddr-xpr:4.0 --nd 1

#hw_emu: host xclbin
#	#add stuff to run

build:  host xclbin


run_system:  build
	./$(host_exe) $(kernel_exe)



