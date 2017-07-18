/*
 * kernel.cl
 *
 *  Created on: Jul 6, 2017
 *      Author: lorenzo.ditucci
 */

#define SIZE 2048

#ifdef __xilinx__
__attribute__((reqd_work_group_size(1,1,1)))
#endif
kernel void vector_add(global uint *a, global uint *b, global uint *out){

uint a_local[SIZE] __attribute((xcl_array_partition(cyclic, 128, 1)));
uint b_local[SIZE] __attribute((xcl_array_partition(cyclic, 128, 1)));
uint out_local[SIZE] __attribute((xcl_array_partition(cyclic, 128, 1)));

#ifdef __xilinx
	__attribute__((xcl_pipelineloop))
#endif
	local_a:for(size_t i = 0; i < SIZE; i++){
		a_local[i] = a[i];
	}

#ifdef __xilinx
	__attribute__((xcl_pipelineloop))
#endif
	local_b:for(size_t i = 0; i < SIZE; i++){
		b_local[i] = b[i];
	}

#ifdef __xilinx
	__attribute__((xcl_pipelineloop))
#endif
	add_loop:for(size_t i = 0; i < SIZE; i++){
		out_local[i] = a_local[i] + b_local[i];
	}

#ifdef __xilinx
	__attribute__((xcl_pipelineloop))
#endif
	global_out:for(size_t i = 0; i < SIZE; i++){
		out[i] = out_local[i];
	}

	return;

}
