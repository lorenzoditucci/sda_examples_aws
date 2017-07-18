#*******************************************************************************
#Vendor: Xilinx 
#Associated Filename: Makefile
#Purpose: Makefile Example for SDAccel
#Author: Jim Wu
#*******************************************************************************
#Copyright (C) 2015 XILINX, Inc.
#
#This file contains confidential and proprietary information of Xilinx, Inc. and 
#is protected under U.S. and international copyright and other intellectual 
#property laws.
#
#DISCLAIMER
#This disclaimer is not a license and does not grant any rights to the materials 
#distributed herewith. Except as otherwise provided in a valid license issued to 
#you by Xilinx, and to the maximum extent permitted by applicable law: 
#(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
#HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
#INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
#FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
#in contract or tort, including negligence, or under any other theory of 
#liability) for any loss or damage of any kind or nature related to, arising under 
#or in connection with these materials, including for any direct, or any indirect, 
#special, incidental, or consequential loss or damage (including loss of data, 
#profits, goodwill, or any type of loss or damage suffered as a result of any 
#action brought by a third party) even if such damage or loss was reasonably 
#foreseeable or Xilinx had been advised of the possibility of the same.
#
#CRITICAL APPLICATIONS
#Xilinx products are not designed or intended to be fail-safe, or for use in any 
#application requiring fail-safe performance, such as life-support or safety 
#devices or systems, Class III medical devices, nuclear facilities, applications 
#related to the deployment of airbags, or any other applications that could lead 
#to death, personal injury, or severe property or environmental damage 
#(individually and collectively, "Critical Applications"). Customer assumes the 
#sole risk and liability of any use of Xilinx products in Critical Applications, 
#subject only to applicable laws and regulations governing limitations on product 
#liability. 
#
#THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
#ALL TIMES.
#
#*******************************************************************************
SHELL = /bin/bash
VPATH = ./

#supported flow: cpu_emu, hw_emu, hw
FLOW=cpu_emu
HOST_SRCS = main.c
HOST_EXE = hello_world
HOST_CFLAGS =  -DFPGA_DEVICE -g -Wall -I${XILINX_OPENCL}/runtime/include/1_2
HOST_LFLAGS = -L${XILINX_OPENCL}/runtime/lib/x86_64 -lxilinxopencl

KERNEL_SRCS = hello_world.cl
KERNEL_DEFS = 
KERNEL_INCS = 
#set target device for XCLBIN
#Set target Xilinx device. Devices below are installed with SDAccel installation:
# xilinx:adm-pcie-7v3:1ddr:2.1 : ADM-PCIE-7V3 Card
# xilinx:adm-pcie-ku3:1ddr:2.1 : ADM-PCIE-KU3 card
# xilinx:adm-pcie-ku3:1ddr:2.1 : ADM-PCIE-KU3 card
XDEVICE=xilinx:adm-pcie-7v3:1ddr:2.1 
KEEP_TEMP=
KERNEL_DEBUG=

CC = g++
CLCC = xocc
CLCC_OPT = --xdevice ${XDEVICE} -o ${XCLBIN} ${KERNEL_DEFS} ${KERNEL_INCS}
ifeq (${FLOW},cpu_emu)
    CLCC_OPT += -t sw_emu
    XCLBIN = hello_world_cpu_emu.xclbin
else ifeq (${FLOW},hw_emu)
    CLCC_OPT += -t hw_emu
    XCLBIN = hello_world_hw_emu.xclbin
else ifeq (${FLOW},hw)
    XCLBIN = hello_world_hw.xclbin
    CLCC_OPT += -t hw
endif

ifeq (${KEEP_TEMP},1)
    CLCC_OPT += -s
endif

ifeq (${KERNEL_DEBUG},1)
    CLCC_OPT += -g
endif

HOST_ARGS = ${XCLBIN}
OBJECTS := $(HOST_SRCS:.c=.o)

.PHONY: all

all: run

run : host xbin
	${HOST_EXE} ${HOST_ARGS}

host: ${HOST_EXE}

xbin : ${XCLBIN}

${HOST_EXE}: ${OBJECTS}
	${CC} ${HOST_LFLAGS} ${OBJECTS} -o $@ 

${XCLBIN}:
	${CLCC} ${CLCC_OPT} ${KERNEL_SRCS}

%.o: %.c
	${CC} ${HOST_CFLAGS} -c $< -o $@

clean:
	${RM} ${HOST_EXE} ${OBJECTS} ${XCLBIN} 


