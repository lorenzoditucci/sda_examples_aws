XOCC=xocc
CC=g++


host_src=./test-cl-2ddr.cpp
host_hdrs=
host_cflags=-D FPGA_DEVICE -g -Wall -I${XILINX_SDX}/runtime/include/1_2 -D C_KERNEL
host_lflags=-L${XILINX_SDX}/runtime/lib/x86_64 -lxilinxopencl

host_exe=host_sw

#kernel
kernel_src=./smithwaterman.cpp
kernel_hdrs=
kernel_flags=
kernel_exe=smithwaterman
kernel_name=smithwaterman

kernel_LDCLFLAGS=--nk smithwaterman:1 \
	--xp param:compiler.preserveHlsOutput=1 \
	--max_memory_ports $(kernel_name) \
	--memory_port_data_width $(kernel_name):512 \
	--xp misc:map_connect=add.kernel.smithwaterman_1.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI \
	--xp misc:map_connect=add.kernel.smithwaterman_1.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI \
	--xp misc:map_connect=add.kernel.smithwaterman_1.M_AXI_GMEM2.core.OCL_REGION_0.M00_AXI 


target_device=xilinx:adm-pcie-ku3:2ddr-xpr:4.0

target=none
ifeq (${target}, sw_emu)
$(info software emulation)
target=sw_emu
else ifeq (${target}, hw_emu)
$(info hardware emulation)
target=hw_emu
else ifeq (${target}, hw)
$(info system build)
target=hw
else
$(info no target selected)
endif

host: $(host_src) $(host_hdrs)
	$(CC) $(host_src) $(host_hdrs) $(host_cflags) $(host_lflags) -o $(host_exe)

xo: 
	$(XOCC) --platform $(target_device) --target $(target) --compile --kernel $(kernel_name) $(kernel_src) $(kernel_LDCLFLAGS)--output $(kernel_exe).xo

xclbin: xo
	$(XOCC) --platform $(target_device) --target $(target) --link --kernel $(kernel_name) $(kernel_exe).xo $(kernel_LDCLFLAGS) --output $(kernel_exe).xclbin

sw_emu: host xclbin
	export XCL_EMULATION_MODE=$(target)
	emconfigutil --xdevice $(target_device) --nd 1
	#emconfigutil -f xilinx:adm-pcie-ku3:2ddr-xpr:4.0 --nd 1

hw_emu: host xclbin
	#add stuff to run

build: host xclbin

run_system: host xclbin
	./$(host_exe) $(kernel_exe)

clean: 
	rm -rf $(host_exe) $(kernel_exe).xo $(kernel_exe).xclbin .Xil emconfig.json 
