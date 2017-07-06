/*
 * main.cpp
 *
 *  Created on: Jul 6, 2017
 *      Author: lorenzo.ditucci
 */

#define SIZE 2048

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

using namespace std;

int
load_file_to_memory(const char *filename, char **result)
{
  int size = 0;
  FILE *f = fopen(filename, "rb");
  if (f == NULL)
  {
    *result = NULL;
    return -1; // -1 means file opening fail
  }
  fseek(f, 0, SEEK_END);
  size = ftell(f);
  fseek(f, 0, SEEK_SET);
  *result = (char *)malloc(size+1);
  if (size != fread(*result, sizeof(char), size, f))
  {
    free(*result);
    return -2; // -2 means file reading fail
  }
  fclose(f);
  (*result)[size] = 0;
  return size;
}


int main (int argc, char *argv[]){

	int err;

	uint *a;
	uint *b;
	uint *out;

	a = new uint[SIZE];
	b = new uint[SIZE];
	out = new uint[SIZE];

	size_t global[2];                   // global domain size for our calculation
	  size_t local[2];                    // local domain size for our calculation

	  cl_platform_id platform_id;         // platform id
	  cl_device_id device_id;             // compute device id
	  cl_context context;                 // compute context
	  cl_command_queue commands;          // compute command queue
	  cl_program program;                 // compute program
	  cl_kernel kernel;                   // compute kernel

	  char cl_platform_vendor[1001];
	  char cl_platform_name[1001];

	  cl_mem input_a;                     // device memory used for the input array
	  cl_mem input_b;                     // device memory used for the input array
	  cl_mem output;                      // device memory used for the output array

	  if (argc != 2){
	    printf("%s <inputfile>\n", argv[0]);
	    return EXIT_FAILURE;
	  }

	for(int i = 0; i < SIZE; i++){
		a[i] = i;
		b[i] = i;
		out[i] = 0;
	}

	// Connect to first platform
	  //
		printf("GET platform \n");
	  err = clGetPlatformIDs(1,&platform_id,NULL);
	  if (err != CL_SUCCESS)
	  {
	    printf("Error: Failed to find an OpenCL platform!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }
		printf("GET platform vendor \n");
	  err = clGetPlatformInfo(platform_id,CL_PLATFORM_VENDOR,1000,(void *)cl_platform_vendor,NULL);
	  if (err != CL_SUCCESS)
	  {
	    printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }
	  printf("CL_PLATFORM_VENDOR %s\n",cl_platform_vendor);
		printf("GET platform name \n");
	  err = clGetPlatformInfo(platform_id,CL_PLATFORM_NAME,1000,(void *)cl_platform_name,NULL);
	  if (err != CL_SUCCESS)
	  {
	    printf("Error: clGetPlatformInfo(CL_PLATFORM_NAME) failed!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }
	  printf("CL_PLATFORM_NAME %s\n",cl_platform_name);

	  // Connect to a compute device
	  //
	  int fpga = 0;
	//#if defined (FPGA_DEVICE)
	  fpga = 1;
	//#endif
		printf("get device \n");
	  err = clGetDeviceIDs(platform_id, fpga ? CL_DEVICE_TYPE_ACCELERATOR : CL_DEVICE_TYPE_CPU,
	                       1, &device_id, NULL);
	  if (err != CL_SUCCESS)
	  {
	    printf("Error: Failed to create a device group!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  // Create a compute context
	  //
		printf("create context \n");
	  context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	  if (!context)
	  {
	    printf("Error: Failed to create a compute context!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  // Create a command commands
	  //
		printf("create queue \n");
	  commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
	  if (!commands)
	  {
	    printf("Error: Failed to create a command commands!\n");
	    printf("Error: code %i\n",err);
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  int status;

	  // Create Program Objects
	  //

	  // Load binary from disk
	  unsigned char *kernelbinary;
	  char *xclbin=argv[1];
	  printf("loading %s\n", xclbin);
	  int n_i = load_file_to_memory(xclbin, (char **) &kernelbinary);
	  if (n_i < 0) {
	    printf("failed to load kernel from xclbin: %s\n", xclbin);
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }
	  size_t n = n_i;
	  // Create the compute program from offline
		printf("create program with binary \n");
	  program = clCreateProgramWithBinary(context, 1, &device_id, &n,
	                                      (const unsigned char **) &kernelbinary, &status, &err);
	  if ((!program) || (err!=CL_SUCCESS)) {
	    printf("Error: Failed to create compute program from binary %d!\n", err);
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  // Build the program executable
	  //
		printf("build program \n");
	  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	  if (err != CL_SUCCESS)
	  {
	    size_t len;
	    char buffer[2048];

	    printf("Error: Failed to build program executable!\n");
	    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
	    printf("%s\n", buffer);
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  // Create the compute kernel in the program we wish to run
	  //
		printf("create kernel \n");
	  kernel = clCreateKernel(program, "vector_add", &err);
	  if (!kernel || err != CL_SUCCESS)
	  {
	    printf("Error: Failed to create compute kernel!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }



		printf("create buffer 0 \n");
	  input_a = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(uint) * SIZE, NULL, NULL);
		printf("create buffer 1 \n");
	  input_b = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(uint) * SIZE, NULL, NULL);
		printf("create buffer 2 \n");
	  output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(uint)* SIZE, NULL, NULL);

	if (!input_a || !input_b || !output)
	  {
	    printf("Error: Failed to allocate device memory!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  // Write our data set into the input array in device memory
	  //
	  //err = clEnqueueWriteBuffer(commands, input_a, CL_TRUE, 0, sizeof(unsigned char) * N , a, 0, NULL, NULL);

		printf("write buffer 0 \n");
	  err = clEnqueueWriteBuffer(commands, input_a, CL_TRUE, 0, sizeof(uint) * SIZE , a, 0, NULL, NULL);
	  if (err != CL_SUCCESS)
	  {
	    printf("Error: Failed to write to source array a!\n");
	    printf("Test failed\n");
	    return EXIT_FAILURE;
	  }

	  err = clEnqueueWriteBuffer(commands, input_b, CL_TRUE, 0, sizeof(uint) * SIZE , b, 0, NULL, NULL);
	  	  if (err != CL_SUCCESS)
	  	  {
	  	    printf("Error: Failed to write to source array a!\n");
	  	    printf("Test failed\n");
	  	    return EXIT_FAILURE;
	  	  }

	  	// Set the arguments to our compute kernel
	  	  //
	  	  err = 0;
	  		printf("set arg 0 \n");
	  	  err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_a);
	  	  if (err != CL_SUCCESS)
	  	  {
	  	    printf("Error: Failed to set kernel arguments 0! %d\n", err);
	  	    printf("Test failed\n");
	  	    return EXIT_FAILURE;
	  	  }
	  		printf("set arg 1 \n");
	  	  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_b);
	  	  if (err != CL_SUCCESS)
	  	  {
	  	    printf("Error: Failed to set kernel arguments 1! %d\n", err);
	  	    printf("Test failed\n");
	  	    return EXIT_FAILURE;
	  	  }
	  		printf("set arg 2 \n");
	  	  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &output);
	  	  if (err != CL_SUCCESS)
	  	  {
	  	    printf("Error: Failed to set kernel arguments 2! %d\n", err);
	  	    printf("Test failed\n");
	  	    return EXIT_FAILURE;
	  	  }

	  	  // Execute the kernel over the entire range of our 1d input data set
	  	  // using the maximum number of work group items for this device
	  	  //
	  	cl_event enqueue_kernel;
//	  	#ifdef C_KERNEL
//	  		printf("LAUNCH task \n");
//	  	  err = clEnqueueTask(commands, kernel, 0, NULL, &enqueue_kernel);
//	  	#else
	  	  global[0] = 1;
	  	  global[1] = 1;
	  	  local[0] = 1;
	  	  local[1] = 1;
	  	  err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
	  	                               (size_t*)&global, (size_t*)&local, 0, NULL, NULL);
//	  	#endif
	  	  if (err)
	  	  {
	  	    printf("Error: Failed to execute kernel! %d\n", err);
	  	    printf("Test failed\n");
	  	    return EXIT_FAILURE;
	  	  }

	  	  // Read back the results from the device to verify the output
	  	  //
	  	  cl_event readevent;
	  		printf("read buffer \n");
	  	  err = clEnqueueReadBuffer( commands, output, CL_TRUE, 0, sizeof(uint) * SIZE, out, 0, NULL, &readevent );
	  	  if (err != CL_SUCCESS)
	  	  {
	  	    printf("Error: Failed to read output array! %d\n", err);
	  	    printf("Test failed\n");
	  	    return EXIT_FAILURE;
	  	  }

	  	  clWaitForEvents(1, &readevent);
	  	  clWaitForEvents(1, &enqueue_kernel);

	  	  for(int i = 0; i < SIZE; i++){
	  		  if(out[i] != i + i){
	  			  printf("Error at index %d! Expected %d, hw: %d \n", i, i+i, out[i]);
	  			  return -1;
	  		  }
	  	  }

	  	  delete [] a;
	  	  delete [] b;
	  	  delete [] out;

	  	  return 0;





}

