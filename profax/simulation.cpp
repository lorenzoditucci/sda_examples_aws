/*
 * Copyright 2016 Lorenzo Di Tucci, Giulia Guidi, Emanuele Del Sozzo
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Freely modified and licensed from http://robotics.stanford.edu/~itayl/mcs/, 
 * with the explicit permission of the autor (Itay Lotan).
*/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
//#include <OpenCL/opencl.h>
#include <CL/opencl.h>
//#include <OpenCL/cl.h>
#include "slist.hpp"
#include "pairtree.hpp"
#include "chain.hpp"
#include "pdb.hpp"
#include "crmsd.hpp"
//#include "cl.hpp"

/*******************************************/
    //  OpenCL

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

    //buffer def
    cl_mem int_array_buff;
    cl_mem real_array_buff;
    cl_mem temp_sum_buff;
    cl_mem type1_m_atypes_dim_buff;
    cl_mem type2_m_atypes_dim_buff;
    cl_mem type1_m_groups_dim_buff;
    cl_mem type2_m_groups_dim_buff;
    cl_mem type1_m_charges_dim_buff;
    cl_mem type2_m_charges_dim_buff;
    cl_mem interaction_buff;
    cl_mem resetTerm_buff;
    cl_mem type1_m_charges_buff;

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




using namespace std;

void generateTorsionMove(vector<ANGLE_CHANGE> & angles, int length, int num_moves,
                         REAL std_angle);
void generateTorsionMove2(vector<ANGLE_CHANGE> & angles, int numres, int num_moves,
                          REAL std_angle, const CChain & chain);
void generateRotamerMove(vector<ROTAMER_CHANGE> & rotamers,
                         int length, int num_rot_changes,
                         const CChain & chain);
bool acceptMove(REAL prev, REAL curr);
float normRand();
REAL attemptMCMove(const vector<ANGLE_CHANGE> & angles,
                   const vector<ROTAMER_CHANGE> & rotamers,
                   CSkiplist & sl, REAL energy, bool bBack, int num_angs);
void storeState(int step, REAL energy, ofstream & fout, CChain & chain);

/*******************************************************************************
 Options:
 
 -I file: The structure input file. The allowed formats is:
 AA name phi-angle psi-angle rotamer index
 The file must have a ".angs" extension.
 -T temp:     The temperature of the simulation in Kelvin.(default is 298)
 -S steps:    The number of trial steps in the MCS. (default is 100,000)
 -N interval: How often to store the current structure. Namely, the number of
 trial steps between save operations (default is 1000).
 -K bbangs:   The number of backbone angles to change simultaneously
 (default is 1)
 -F numrot    How many rotamer step to make per backbone step (default is 5)
 -R rots:     The number of rotamers to change each time (default is 1)
 -O output:   The name of the output file prefix (including directory name).
 The default is to only save final results.
 -U rseed:    The random seed to be used.
 -A angle:    The std of each angular change.
 -C rfile     Use this structure as a reference for all cRMSD computations
 -H dir:      The directory where the executable resides.
 -X tempdir   The directory where intermediate results will be saved.
 
 *******************************************************************************/

char ifile[100] = "";
char rfile[100] = "";
REAL temper = T_0;
REAL temperature;
int numSteps = 1000;
int interval = 1000;
int bbangs = 1;
int numrot = 5;
int rotangs = 1;
char ofile[100] = "";
int ifType;
int rseed = 0;
bool noOutput = false;
REAL stdAngle = 10*M_PI/180.0;
bool crmsd;
bool refStruct = false;
char dir[100] = ".";
char tempdir[100] = "";
REAL angle_factor = 1.0;

POSITIONS pos0, posCurr;;

void parseNextItem(int argc, char ** argv, int & i)
{
    if (strncmp(argv[i], "-I", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No input file name provided!" << endl;
            exit(0);
        }
        
        strcpy(ifile, argv[++i]);
        
    }
    else if (strncmp(argv[i], "-X", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No temporary directory name provided!" << endl;
            exit(0);
        }
        
        strcpy(tempdir, argv[++i]);
        
    }
    else if (strncmp(argv[i], "-H", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No home directory was provided!" << endl;
            exit(0);
        }
        
        strcpy(dir, argv[++i]);
        
    }
    else if (strncmp(argv[i], "-C", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "A pdb files need to be provided for -C option!" << endl;
            exit(0);
        }
        
        refStruct = true;
        strcpy(rfile, argv[++i]);
    }
    else if (strncmp(argv[i], "-T", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No system temperature provided!" << endl;
            exit(0);
        }
        temper = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "-S", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No  number of trial steps provided!" << endl;
            exit(0);
        }
        numSteps = atoi(argv[++i]);
    }
    
    else if (strncmp(argv[i], "-N", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No interval size provided!" << endl;
            exit(0);
        }
        interval = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-K", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No number of backbone angles to change per step provided!" << endl;
            exit(0);
        }
        bbangs = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-F", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "The number of rotamer moves to make per backbone move!" << endl;
            exit(0);
        }
        numrot = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-R", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No number of rotamers to change per rotamer step provided!" << endl;
            exit(0);
        }
        rotangs = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-O", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No output file prefix provided!" << endl;
            exit(0);
        }
        strcpy(ofile, argv[++i]);
    }
    else if (strncmp(argv[i], "-U", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No random seed provided!" << endl;
            exit(0);
        }
        rseed = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-A", 2) == 0)
    {
        if (argc < i + 2)
        {
            cout << "No std of backbone angular change provided!" << endl;
            exit(0);
        }
        stdAngle = atof(argv[++i])*M_PI/180.0;
    }
    else
    {
        cout << "Bad command line input!!!" << endl;
        exit(0);
    }
    
    i++;
}

void parseCommandLine(int argc, char ** argv)
{
    int i = 1;
    while (i < argc)
        parseNextItem(argc, argv, i);
    
    char * ext = strrchr(ifile, '.');
    if (!ext)
    {
        cout << "No input file extension! " << ifile << " " << endl;
        exit(0);
    }
    
    if (strcmp(".angs", ext) == 0)
        ifType = ANGS_FILE;
    else
    {
        cout << "Bad input file extension, wrong type of file!" << endl;
        exit(0);
    }
    
    if (refStruct)
    {
        char * ext1 = strrchr(rfile, '.');
        
        if (strcmp(".pdb", ext1) != 0)
        {
            cout << "Input refernce structure file needs to be a PDB file" << endl;
            exit(0);
        }
    }
    
    if (rotangs == 0)
        numrot = 0;
    
    if (strcmp(ofile, "") == 0)
        noOutput = true;
    
    temperature = 1.0/(Kb * temper);
}

void printOptions(ostream & out)
{
    out << "The input file is: " << ifile << endl;
    out << "The simulation will have " << numSteps << " trial steps" << endl;
    out << "Each step consists of a backbone move followed by " << numrot
    << " sidechain moves" << endl;
    out << "Each backbone move consists of changing " << bbangs
    << " backbone angles with a StDev of " << stdAngle*180.0/M_PI << endl;
    if (numrot > 0)
        out << "Each rotamer move consists of changing " << rotangs << " rotamers"
        << " selected at random." << endl;
    out << "The temperature of the simulation is " << temper << endl;
    out << "Intermediate results will be saved every " << interval
    << " trial backbone steps" << endl;
    if (!noOutput)
        out << "Output will be written to " << ofile << " with suitable extensions" << endl;
    else
        out << "Only final output will be saved in the current directory" << endl;
    out << "The random seed is " << rseed << endl;
    if (refStruct)
        out << "The reference structure is: " << rfile << endl;
}

int main(int argc, char ** argv)
{


    parseCommandLine(argc, argv);
    
    printOptions(cout);
    Initialize(dir);
    
    ofstream fout;
    srand48(rseed);
    
    // Open main output file
    char buf[100];
    sprintf(buf, "%s/%s.out", dir, ofile);
    fout.open(buf);
    
    if (!fout.is_open())
    {
        cout << "Could not open output file " << buf << endl;
        exit(1);
    }
    
    // Store the command lline options in a file.
    sprintf(buf, "%s/%s.opt", dir, ofile);
    ofstream fopt(buf);
    printOptions(fopt);
    
    // Create the chain from the description in the file.
    char f[100];
    sprintf(f, "%s/%s", dir, ifile);
    CChain chain(f, ifType);
    
    // create chaintree
    CSkiplist sl(chain);
    
    vector<ANGLE_CHANGE> angles(bbangs + 1);
    vector<ROTAMER_CHANGE> rotamers(rotangs + 1);
    
    vector<ANGLE_CHANGE> no_angle_change(1);
    no_angle_change[0].m_index = chain.getLength() + 1;
    vector<ROTAMER_CHANGE> no_rotamer_change(1);
    no_rotamer_change[0].m_index = chain.getLength() + 1;
   


 int i;

    

    i = 0;
    
    int err = clGetPlatformIDs(1,&platform_id,NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: Failed to find an OpenCL platform!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  err = clGetPlatformInfo(platform_id,CL_PLATFORM_VENDOR,1000,(void *)cl_platform_vendor,NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  printf("CL_PLATFORM_VENDOR %s\n",cl_platform_vendor);
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
#if defined (FPGA_DEVICE)
  fpga = 1;
#endif

  printf("FPGA is %d \n", fpga);
  err = clGetDeviceIDs(platform_id, fpga ? CL_DEVICE_TYPE_ACCELERATOR : CL_DEVICE_TYPE_CPU,
                       1, &device_id, NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: Failed to create a device group! %d \n", err);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  
  // Create a compute context 
  //
  context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
  if (!context)
  {
    printf("Error: Failed to create a compute context!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }

  // Create a command commands
  //
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
  char *xclbin="computePairEnergy.xclbin";
   // char *xclbin="kernel.cl";
  printf("loading %s\n", xclbin);
  int n_i = load_file_to_memory(xclbin, (char **) &kernelbinary);
  if (n_i < 0) {
    printf("failed to load kernel from xclbin: %s\n", xclbin);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  size_t n = n_i;
  // Create the compute program from offline
  program = clCreateProgramWithBinary(context, 1, &device_id, &n,
                                      (const unsigned char **) &kernelbinary, &status, &err);
    //program = clCreateProgramWithSource(context, 1, (const char **) & kernelbinary, NULL, &err);
  if ((!program) || (err!=CL_SUCCESS)) {
    printf("Error: Failed to create compute program from binary %d!\n", err);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }

  // Build the program executable
  //
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
  kernel = clCreateKernel(program, "computePairEnergy_k", &err);
  if (!kernel || err != CL_SUCCESS)
  {
    printf("Error: Failed to create compute kernel!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }

  //create the buffers..
  cout << "creating buffers..." << endl;
  int_array_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int) *(8 + MAX_TYPE1_M_ATYPES_DIM + MAX_TYPE2_M_ATYPES_DIM + MAX_TYPE1_M_GROUPS_DIM + MAX_TYPE2_M_GROUPS_DIM), NULL, NULL);
  real_array_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL) * ( 9 + 3 + 3 + 3 + 1 + 1 + 9 + 12*3 + 3 + 9 + 12*3 + 3 + 12*3 + 12*3 + MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE + MAX_TYPE1_M_CHARGES_DIM + MAX_TYPE2_M_CHARGES_DIM), NULL, NULL);
  temp_sum_buff = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(REAL), NULL, NULL);
 
	cout << "random print" << endl; 
  type1_m_atypes_dim_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
  type2_m_atypes_dim_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
  
  type1_m_groups_dim_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
  type2_m_groups_dim_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
  
  type1_m_charges_dim_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
  type2_m_charges_dim_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);

  interaction_buff = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(int), NULL, NULL);  
  resetTerm_buff = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(short), NULL, NULL);  

  type1_m_charges_buff = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(REAL) * MAX_TYPE1_M_CHARGES_DIM, NULL, NULL);

  cout << "buffers created " << endl;



        
    // Compute the initial energy of the structure
    REAL energy = sl.computeEnergy(no_angle_change);
    cout << "The energy of the initial conformation is " << energy << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    
    // Store the initial conformations.
    char cc[100];
    if (!noOutput)
    {
        sprintf(cc, "%s/%s_%06d.pdb", tempdir, ofile, 0);
        chain.storeCoordinates(cc);
        sprintf(cc, "%s/%s_%06d.angs", tempdir, ofile, 0);
        chain.storeAngsStyle(cc);
    }
    
    // Load the reference structure.
    if (refStruct)
    {
        sprintf(f, "%s/%s", dir, rfile);
        loadCas(f, pos0);
    }
    else
        // If none given use the initial structure as reference.
        chain.getCaPositions(pos0);
    
    // The main loop of the simulation.
    //**************************************************************
    for (i = 1; i <= numSteps; i++)
    {
        // Perform the backbone move.
        generateTorsionMove(angles, chain.getLength(), bbangs, stdAngle);
        
        energy = attemptMCMove(angles, no_rotamer_change, sl, energy, true,
                               bbangs);
        
        // Perform all the rotamer moves
        for (int j = 0; j < numrot; j++)
        {
            generateRotamerMove(rotamers, chain.getLength(), rotangs, chain);
            energy = attemptMCMove(no_angle_change, rotamers, sl, energy,
                                   false, -1);
        }
        
        // Store the current state every 'interval' number of steps.
        if (i % interval == 0)
            storeState(i, energy, fout, chain);
    }
    //********************************************************************
    
    if (noOutput)
    {
        chain.storeCoordinates("output.pdb");
        chain.storeAngsStyle("output.angs");
    }

    //clean!

    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    clReleaseMemObject(int_array_buff);
    clReleaseMemObject(real_array_buff);
    clReleaseMemObject(temp_sum_buff);
    clReleaseMemObject(type1_m_atypes_dim_buff);
    clReleaseMemObject(type2_m_atypes_dim_buff);
    clReleaseMemObject(type1_m_groups_dim_buff);
    clReleaseMemObject(type2_m_groups_dim_buff);
    clReleaseMemObject(type1_m_charges_dim_buff);
    clReleaseMemObject(type2_m_charges_dim_buff);
    clReleaseMemObject(interaction_buff);
    clReleaseMemObject(resetTerm_buff);
    
    
    return 0;
}

REAL attemptMCMove(const vector<ANGLE_CHANGE> & angles,
                   const vector<ROTAMER_CHANGE> & rotamers,
                   CSkiplist & sl, REAL energy, bool bBack, int num_angs)
{
    sl.makeMove(angles, rotamers);
    REAL newEnergy;
    
    if (sl.findSelfClash())
    {
        sl.undoLastMove();
        return energy;
    }
    else
    {
        newEnergy = sl.computeEnergy(angles);
        
        if (!acceptMove(energy, newEnergy))
        {
            sl.undoLastMove();
            return energy;
        }
        else
            return newEnergy;
    }
}

// The acceptance rule for the MC simulation.
bool acceptMove(REAL prev, REAL curr)
{
    REAL diff = prev - curr;
    
    if (diff >= 0)
        return true;
    else if (diff*temperature < -15)
        return false;
    else
    {
        REAL e = exp(diff*temperature);
        REAL r = drand48();
        if (r < e)
            return true;
        else
            return false;
    }
}

// A random number generator that picks numbers from a normal distribution.
float normRand()
{
    static float saved;
    static bool bSaved = false;
    
    float v1, v2, rsq;
    if (!bSaved)
    {
        do {
            v1 = 2.0 * drand48() - 1.0;
            v2 = 2.0 * drand48() - 1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        float fac = sqrt(-2.0*log(rsq)/rsq);
        saved = v1*fac;
        bSaved = true;
        return fac*v2;
    }
    else
    {
        bSaved = false;
        return saved;
    }
}


void generateTorsionMove(vector<ANGLE_CHANGE> & angles, int length, int num_moves,
                         REAL std_angle)
{
    int k;
    int used[num_moves];
    
    // Generate torsion moves
    for (k = 0; k < num_moves; k++)
    {
        bool found = false;
        
        // The first phi angle and last psi angle are not changed.
        int index = lrand48() % (length - 3) + 1;
        
        for (int p = 0; p < k; p++)
        {
            if (used[p] == index)
            {
                found = true;
                break;
            }
        }
        if (found)
        {
            k--;
            continue;
        }
        
        used[k] = index;
        angles[k].m_index = index;
        angles[k].m_angle = normRand() * std_angle;
    }
    
    // Insert stopper.
    angles[k].m_index = length + 1;
    
    sort(angles.begin(), angles.end(), ANGLE_CHANGE_COMP());
}

void generateRotamerMove(vector<ROTAMER_CHANGE> & rotamers,
                         int length, int num_rot_changes,
                         const CChain & chain)
{
    int k;
    
    // Generate rotamer moves.
    for (k = 0; k < num_rot_changes; k++)
    {
        bool found = false;
        
        // Only odd numbered links contain side-chains
        int index = (lrand48() % ((length - 1)/2))*2 + 1;
        
        for (int p = 0; p < k; p++)
            if (rotamers[p].m_index == index)
            {
                found = true;
                break;
            }
        
        if (found)
        { 
            k--;
            continue;
        }
        
        int n = SIDECHAIN::m_aalist[chain.getLink(index)->getType()]->m_nRotamers;
        
        int ri;
        if (n > 1)
        {
            ri = (lrand48() % (n - 1));
            if (ri >= chain.getLink(index)->getRotIndex())
                ri++;
        }
        else
        {
            k--;
            continue;
        }
        
        rotamers[k].m_index = index;     
        rotamers[k].m_rotIndex = ri;
    }
    
    // Insert stopper
    rotamers[k].m_index = length + 1;
    
    sort(rotamers.begin(), rotamers.end(), ROTAMER_CHANGE_COMP());  
}

// Store the current state of the simulation.
void storeState(int step, REAL energy, ofstream & fout, CChain & chain)
{
    char cc[100];
    
    if (!noOutput)
    {
        int index = step/interval;
        
        sprintf(cc, "%s/%s_%06d.pdb", tempdir, ofile, index);
        chain.storeCoordinates(cc);
        sprintf(cc, "%s/%s_%06d.angs", tempdir, ofile, index);
        chain.storeAngsStyle(cc);
    }
    
    chain.getCaPositions(posCurr);
    REAL rot[3][3], trans[3];
    REAL rmsd = CRMSD(pos0, posCurr, rot, trans);
    
    fout << energy << " " << rmsd <<  endl;
}


