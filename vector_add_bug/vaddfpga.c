/*
 * vaddfpga.c

 *
 *  Created on: Jul 11, 2017
 *      Author: emanuele.delsozzo
 */

#include "vadd.h"

void krnl_vadd(int a[LENGTH], int b[LENGTH], int c[LENGTH]){
#pragma HLS INTERFACE m_axi port=a offset=slave
#pragma HLS INTERFACE m_axi port=b offset=slave
#pragma HLS INTERFACE m_axi port=c offset=slave

#pragma HLS INTERFACE s_axilite port=a bundle=control
#pragma HLS INTERFACE s_axilite port=b bundle=control
#pragma HLS INTERFACE s_axilite port=c bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control


	for(int i = 0; i < LENGTH; i++){
	    c[i] = a[i] + b[i] + VAL;
	  }



}


