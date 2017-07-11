# SHA1 Application
COMMON_REPO =./

include $(COMMON_REPO)/utility/boards.mk
include $(COMMON_REPO)/libs/logger/logger.mk
include $(COMMON_REPO)/libs/cmdparser/cmdparser.mk
include $(COMMON_REPO)/libs/opencl/opencl.mk

# hello Host Application
smithwaterman_SRCS=./test-cl-2ddr.cpp  $(logger_SRCS) $(cmdparser_SRCS) $(xcl_SRCS)
smithwaterman_HDRS= $(logger_SRCS) $(cmdparser_HDRS) $(xcl_HDRS)
smithwaterman_CXXFLAGS=-std=gnu++0x -I./ $(opencl_CXXFLAGS) $(logger_CXXFLAGS) $(cmdparser_CXXFLAGS) $(xcl_CXXFLAGS) -O3 -Wall -D FPGA_DEVICE -D C_KERNEL -I$(XILINX_SDACCEL)/Vivado_HLS/include
smithwaterman_LDFLAGS=$(opencl_LDFLAGS) -lrt

EXES=smithwaterman

# SHA1 Kernel
krnl_smithwaterman_SRCS=./smithwaterman.cpp
krnl_smithwaterman_CLFLAGS=-I./ -k smithwaterman
krnl_smithwaterman_NDEVICES=xilinx:xil-accel-rd-vu9p xilinx:xil-accel-rd-ku115:4ddr-xpr xilinx:adm-pcie-ku3:2ddr-xpr xilinx:adm-pcie-7v3:1ddr

XOS=krnl_smithwaterman

# SHA1 xclbin
krnl_smithwaterman_XOS=krnl_smithwaterman
krnl_smithwaterman_LDCLFLAGS=--nk smithwaterman:1 \
	--xp param:compiler.preserveHlsOutput=1 \
	--max_memory_ports \
	--memory_port_data_width [$(XOS)]:512 \
	--xp misc:map_connect=add.kernel.smithwaterman.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI \
	--xp misc:map_connect=add.kernel.smithwaterman.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI \
	--xp misc:map_connect=add.kernel.smithwaterman.M_AXI_GMEM2.core.OCL_REGION_0.M00_AXI 

krnl_smithwaterman_NDEVICES=xilinx:xil-accel-rd-vu9p:4ddr-xpr xilinx:xil-accel-rd-ku115:4ddr-xpr xilinx:adm-pcie-ku3:2ddr-xpr xilinx:adm-pcie-7v3:1ddr


XCLBINS=krnl_smithwaterman

# check
check_EXE=smithwaterman
check_XCLBINS=krnl_smithwaterman
check_NDEVICES=$(krnl_smithwaterman_NDEVICES)

# Pattern to allow cmdline to determine xclbin name
define mk_args
#check_$(1)_$(call sanitize_dsa,$(2))_ARGS = -t 20 -p Xilinx -d $(2) -k ./xclbin/krnl_smithwaterman.$(1).$(call sanitize_dsa,$(2)).xclbin
check_$(1)_$(call sanitize_dsa,$(2))_ARGS = ./xclbin/krnl_smithwaterman.$(1).$(call sanitize_dsa,$(2)).xclbin
endef

#DEVICES=xilinx:adm-pcie-ku3:2ddr-xpr:4.0
DEVICES=xilinx:adm-pcie-ku3:2ddr-xpr:4.0


$(foreach target,$(TARGETS),$(foreach device,$(DEVICES),$(eval $(call mk_args,$(target),$(device)))))

CHECKS=check

include $(COMMON_REPO)/utility/rules.mk

