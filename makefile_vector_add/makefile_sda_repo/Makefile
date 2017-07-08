# SHA1 Application
COMMON_REPO =./

include $(COMMON_REPO)/utility/boards.mk
include $(COMMON_REPO)/libs/logger/logger.mk
include $(COMMON_REPO)/libs/cmdparser/cmdparser.mk
include $(COMMON_REPO)/libs/opencl/opencl.mk

# hello Host Application
vec_add_SRCS=main.cpp  $(logger_SRCS) $(cmdparser_SRCS) $(xcl_SRCS)
vec_add_HDRS= $(logger_SRCS) $(cmdparser_HDRS) $(xcl_HDRS)
vec_add_CXXFLAGS=-std=gnu++0x -I./src/ $(opencl_CXXFLAGS) $(logger_CXXFLAGS) $(cmdparser_CXXFLAGS) $(xcl_CXXFLAGS)
vec_add_LDFLAGS=$(opencl_LDFLAGS) -lrt

EXES=vec_add

# SHA1 Kernel
krnl_vec_add_SRCS=./kernel.cl
krnl_vec_add_HDRS=
krnl_vec_add_CLFLAGS=-I./src/
krnl_vec_add_NDEVICES=xilinx:xil-accel-rd-ku115:4ddr-xpr:3.3 xilinx:adm-pcie-ku3:2ddr-xpr:4.0

XOS=krnl_vec_add

# SHA1 xclbin
krnl_vec_add_XOS=krnl_vec_add
krnl_vec_add_LDCLFLAGS=--nk dev_vec_add_update:2
krnl_vec_add_NDEVICES=xilinx:xil-accel-rd-ku115:4ddr-xpr:3.3 xilinx:adm-pcie-ku3:2ddr-xpr:4.0

XCLBINS=krnl_vec_add

# check
check_EXE=vec_add
check_XCLBINS=krnl_vec_add
check_NDEVICES=$(krnl_vec_add_NDEVICES)

# Pattern to allow cmdline to determine xclbin name
define mk_args
check_$(1)_$(call sanitize_dsa,$(2))_ARGS = -t 20 -p Xilinx -d $(2) -k ./xclbin/krnl_vec_add.$(1).$(call sanitize_dsa,$(2)).xclbin
endef

DEVICES=xilinx:adm-pcie-ku3:2ddr-xpr:4.0

$(foreach target,$(TARGETS),$(foreach device,$(DEVICES),$(eval $(call mk_args,$(target),$(device)))))

CHECKS=check

include $(COMMON_REPO)/utility/rules.mk

