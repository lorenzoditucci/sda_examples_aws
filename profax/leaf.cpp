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


#include <iostream>

#include "leaf.hpp"
#include "pairtree.hpp"
#include "slist.hpp"
#include "leafCL.hpp"
#include "spheres.hpp"
#include "kernel.h"
#include "RectDistCL.h"
#include <CL/opencl.h>

/*******************************************/
//  OpenCL

extern size_t global[2];                   // global domain size for our calculation
extern size_t local[2];                    // local domain size for our calculation

extern cl_platform_id platform_id;         // platform id
extern cl_device_id device_id;             // compute device id
extern cl_context context;                 // compute context
extern cl_command_queue commands;          // compute command queue
extern cl_program program;                 // compute program
extern cl_kernel kernel;                   // compute kernel

extern char cl_platform_vendor[1001];
extern char cl_platform_name[1001];

//buffer def
extern cl_mem int_array_buff;
extern cl_mem real_array_buff;
extern cl_mem temp_sum_buff;
extern cl_mem type1_m_atypes_dim_buff;
extern cl_mem type2_m_atypes_dim_buff;
extern cl_mem type1_m_groups_dim_buff;
extern cl_mem type2_m_groups_dim_buff;
extern cl_mem type1_m_charges_dim_buff;
extern cl_mem type2_m_charges_dim_buff;
extern cl_mem interaction_buff;
extern cl_mem resetTerm_buff;
extern cl_mem type1_m_charges_buff;


short resetTerm[1];
int interaction[1];
REAL tempSum[1];


REAL CLeaf::m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE];

// Switch the rotamer coordinates held by this leaf to different
// precomputed values.
void CLeaf::changeRotamer(int rotIndex)
{
    assert(rotIndex < SIDECHAIN::m_aalist[getType()]->m_nRotamers);
    
    m_undoRotIndex = m_rotIndex;
    m_rotIndex = rotIndex;
    
    const ROTAMER & rot =
    SIDECHAIN::m_aalist[getType()]->getRotamer(rotIndex, m_rotType);
    
    // Set the new coordinates as well as their BV and energy, which were
    // precomputed.
    m_positions = rot.m_positions;
    m_bv = rot.m_bv;
    m_energy = rot.m_energy;
}

// Compute all pairs of distances between the atoms of this leaf and
// the given leaf.
void CLeaf::computeDistances(CLeaf * pLeaf, const REAL rot[3][3],
                             const REAL trans[3])
{
    int size1 = SIDECHAIN::m_aalist[getType()]->m_size;
    int size2 = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
    
    REAL cen[3], dist[3], vec1[3], vec2[3];
    
    // If this is a PRO backbone, we need to add the Cd atom because
    // it is part of a group with the Ca atom for Electrostatic purposes
    if (getType() == BBP)
    {
        assert(getNext());
        assert(getNext()->getType() == PRO);
        MxVpV(vec1, m_rotate, getNext()->getPositions()[2], m_translate);
    }
    
    // If this is a PRO backbone, we need to add the Cd atom because
    // it is part of a group with the Ca atom for Electrostatic purposes
    if (pLeaf->getType() == BBP)
    {
        assert(pLeaf->getNext());
        assert(pLeaf->getNext()->getType() == PRO);
        
        REAL temp[3];
        MxVpV(temp, pLeaf->m_rotate, pLeaf->getNext()->getPositions()[2],
              pLeaf->m_translate);
        MxVpV(vec2, rot, temp, trans);
    }
    
    // Compute te distances between all pairs of atoms.
    for(int j = 0; j < size2; j++)
    {
        MxVpV(cen, rot, pLeaf->getPositions()[j], trans);
        
        for (int i = 0; i < size1; i++)
        {
            VmV(dist, cen, getPositions()[i]);
            m_distances[i][j] = Vlength2(dist);
            
            // Add distances to the Cd of the second node (if type is BBP)
            if (pLeaf->getType() == BBP)
            {
                VmV(dist, vec2, getPositions()[i]);
                m_distances[i][size2] = Vlength2(dist);
            }
        }
        
        // Add distances to the Cd of the first node (if type is BBP)
        if (getType() == BBP)
        {
            VmV(dist, cen, vec1);
            m_distances[size1][j] = Vlength2(dist);
        }
    }
    
    // Add distance between the Cd of the first and second nodes (both BBPs)
    if (getType() == BBP && pLeaf->getType() == BBP)
    {
        VmV(dist, vec2, vec1);
        m_distances[size1][size2] = Vlength2(dist);
    }
}


void CLeaf::computePairEnergy(CNode * pNode, const REAL rot[3][3],
                              const REAL trans[3], CTerm * term,
                              bool bSeparated)
{
    CLeaf * pLeaf = (CLeaf*) pNode;
    const CSphere * s1 = (CSphere*) this->getBV();
    const CSphere * s2 = (CSphere*) pLeaf->getBV();
    
    COORDS pLeaf_next_pos;
    
    if (pLeaf->getNext()) {
        pLeaf_next_pos = (COORDS)pLeaf->getNext()->getPositions();
    }else{
        pLeaf_next_pos = (COORDS)pLeaf->getPositions();
    }
    
    
    //std::cout << pLeaf->getNext() << "\n";
    
    //instead of this one, need to call the kernel one!
    //computePairEnergyCL(getType(), pLeaf->getType(), pLeaf->isLeaf(), rot, trans, bSeparated, s1->m_center, s1->m_center, s1->m_rad, s2->m_rad, SIDECHAIN::m_aalist[getType()]->m_size, SIDECHAIN::m_aalist[pLeaf->getType()]->m_size, m_rotate, getNext()->getPositions(), m_translate, pLeaf->m_rotate, pLeaf_next_pos, pLeaf->m_translate, pLeaf->getPositions(), getPositions(), m_distances, getIndex(), pLeaf->getIndex(), SIDECHAIN::m_aalist[getType()]->m_aTypes, SIDECHAIN::m_aalist[pLeaf->getType()]->m_aTypes, SIDECHAIN::m_aalist[getType()]->m_nGroups, SIDECHAIN::m_aalist[pLeaf->getType()]->m_nGroups, SIDECHAIN::m_aalist[getType()]->m_groups, SIDECHAIN::m_aalist[pLeaf->getType()]->m_groups, SIDECHAIN::m_aalist[getType()]->m_charges, SIDECHAIN::m_aalist[pLeaf->getType()]->m_charges);
    
    //input arrays definition
    static int * int_array = NULL;
    static REAL * real_array = NULL;
    if(int_array == NULL)
        int_array = (int *)malloc(sizeof(int) *(8 + MAX_TYPE1_M_ATYPES_DIM + MAX_TYPE2_M_ATYPES_DIM + MAX_TYPE1_M_GROUPS_DIM + MAX_TYPE2_M_GROUPS_DIM) );
    if(real_array == NULL)
        real_array = (REAL *)malloc(sizeof(REAL) * ( 9 + 3 + 3 + 3 + 1 + 1 + 9 + 12*3 + 3 + 9 + 12*3 + 3 + 12*3 + 12*3 + MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE + MAX_TYPE1_M_CHARGES_DIM + MAX_TYPE2_M_CHARGES_DIM));
    
    
    //fill up the arrays
    int_array[0] = getType();
    int_array[1] = pLeaf->getType();
    int_array[2] = SIDECHAIN::m_aalist[getType()]->m_size;
    int_array[3] = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
    int_array[4] = this->getIndex();
    int_array[5] = pLeaf->getIndex();
    
    int index = 6;
    for(int i = index, s_i = 0; i < index + SIDECHAIN::m_aalist[this->getType()]->m_size; i++, s_i++) int_array[i] = SIDECHAIN::m_aalist[this->getType()]->m_aTypes[s_i];
    index += SIDECHAIN::m_aalist[this->getType()]->m_size;
    for(int i = index, s_i = 0; i < index + SIDECHAIN::m_aalist[pLeaf->getType()]->m_size; i++, s_i++) int_array[i] = SIDECHAIN::m_aalist[pLeaf->getType()]->m_aTypes[s_i];
    index += SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
    int_array[index] = SIDECHAIN::m_aalist[this->getType()]->m_nGroups;
    index++;
    int_array[index] =  SIDECHAIN::m_aalist[pLeaf->getType()]->m_nGroups;
    index++;
    for(int i = index, s_i = 0; i < index + SIDECHAIN::m_aalist[this->getType()]->m_nGroups; i++, s_i++) int_array[i] = SIDECHAIN::m_aalist[this->getType()]->m_groups[s_i];
    index += SIDECHAIN::m_aalist[this->getType()]->m_nGroups;
    for(int i = index, s_i = 0; i < index + SIDECHAIN::m_aalist[pLeaf->getType()]->m_nGroups; i++, s_i++) int_array[i] = SIDECHAIN::m_aalist[pLeaf->getType()]->m_groups[s_i];
    index += SIDECHAIN::m_aalist[pLeaf->getType()]->m_nGroups;
    
    if(index >= 8 + MAX_TYPE1_M_ATYPES_DIM + MAX_TYPE2_M_ATYPES_DIM + MAX_TYPE1_M_GROUPS_DIM + MAX_TYPE2_M_GROUPS_DIM){
        cout << "ERROR: out of bound on int_array" << endl;;
    }
    
    //buffer REAL
    index = 0;
    //for(int i = index, s_i = 0; i < index + 9; i++, s_i++) real_array[i] = rot[s_i];
    int t_i = 0;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            real_array[t_i] = rot[i][j];
            t_i ++;
        }
    }
    index += 9;
    for(int i = index, s_i = 0; i < index + 3; i++, s_i++) real_array[i] = trans[s_i];
    index += 3;
    for(int i = index, s_i = 0; i < index + 3; i++, s_i++) real_array[i] = s1->m_center[s_i];
    index +=3;
    for(int i = index, s_i = 0; i < index + 3; i++, s_i++) real_array[i] = s2->m_center[s_i];
    index +=3;
    real_array[index] = s1->m_rad;
    index ++;
    real_array[index] = s2->m_rad;
    index ++;
    //for(int i = index, s_i = 0; i < index + 9; i++, s_i++) real_array[i] = this->m_rota[s_i];
    for(int i = index, s_i = 0; s_i < 3; s_i++){
        for(int s_j = 0; s_j < 3; s_j++){
            real_array[i] = this->m_rotate[s_i][s_j];
            i++;
        }
        
    }
    index += 9;
    //not sure of the real dimensions of this one....
    for(int i = index, s_i = 0; s_i < 12; s_i++){
        for(int s_j = 0; s_j < 3; s_j++){
            real_array[i] = this->getNext()->getPositions()[s_i][s_j];
            i++;
        }
    }
    //for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++) real_array[i] = CLeaf_next_m_positions[s_i];
    index += 12 * 3;
    for(int i = index, s_i = 0; i < index + 3; i++, s_i++) real_array[i] = this->m_translate[s_i];
    index += 3;
    
    for(int i = index, s_i = 0; s_i < 3; s_i++){
        for(int s_j = 0; s_j < 3; s_j++){
            real_array[i] = pLeaf->m_rotate[s_i][s_j];
            i++;
        }
        
    }
    //for(int i = index, s_i = 0; i < index + 9; i++, s_i++) real_array[i] = pNode_m_rotate[s_i];
    index += 9;
    
    //not sure of the real dimensions of this one....
    if(pLeaf->getNext() != NULL){
        for(int i = index, s_i = 0; s_i < 12; s_i++){
            for(int s_j = 0; s_j < 3; s_j++){
                //std::cout << "indici " << s_i << " - " << s_j << std::endl;
                real_array[i] = pLeaf->getNext()->getPositions()[s_i][s_j];
                i++;
            }
        }
    }else{
        for(int i = index, s_i = 0; s_i < 12; s_i++){
            for(int s_j = 0; s_j < 3; s_j++){
                //std::cout << "indici " << s_i << " - " << s_j << std::endl;
                real_array[i] = pLeaf->getPositions()[s_i][s_j];
                i++;
            }
        }
    }
    
    //for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++) real_array[i] = pNode_next_m_positions[s_i];
    index += 12 * 3;
    for(int i = index, s_i = 0; i < index + 3; i++, s_i++) real_array[i] = pLeaf->m_translate[s_i];
    index += 3;
    
    //not sure of the real dimensions of this one....
    for(int i = index, s_i = 0; s_i < 12; s_i++){
        for(int s_j = 0; s_j < 3; s_j++){
            real_array[i] = pLeaf->getPositions()[s_i][s_j];
            i++;
        }
    }
    
    //for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++) real_array[i] = pNode_m_positions[s_i];
    index += 12 * 3;
    
    //not sure of the real dimensions of this one....
    for(int i = index, s_i = 0; s_i < 12; s_i++){
        for(int s_j = 0; s_j < 3; s_j++){
            real_array[i] = this->getPositions()[s_i][s_j];
            i++;
        }
    }
    //for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++) real_array[i] = CLeaf_m_positions[s_i];
    index += 12 * 3;
    
    
    for(int i = index, s_i = 0; s_i < MAX_ROTAMER_SIZE; s_i++){
        for(int s_j = 0; s_j < MAX_ROTAMER_SIZE; s_j++){
            real_array[i] = this->m_distances[s_i][s_j];
            i++;
        }
    }
    //for(int i = index, s_i = 0; i < index + MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE; i++, s_i++) real_array[i] = CLeaf_m_distances[s_i];
    index += MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE;
    for(int i = index, s_i = 0; i < index + SIDECHAIN::m_aalist[getType()]->m_size; i++, s_i++) real_array[i] = SIDECHAIN::m_aalist[this->getType()]->m_charges[s_i];
    index += SIDECHAIN::m_aalist[getType()]->m_size;
    for(int i = index, s_i = 0; i < index + SIDECHAIN::m_aalist[pLeaf->getType()]->m_size; i++, s_i++) real_array[i] = SIDECHAIN::m_aalist[pLeaf->getType()]->m_charges[s_i];
    index += SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
    
    if(index >= 9 + 3 + 3 + 3 + 1 + 1 + 9 + 12*3 + 3 + 9 + 12*3 + 3 + 12*3 + 12*3 + MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE + MAX_TYPE1_M_CHARGES_DIM + MAX_TYPE2_M_CHARGES_DIM){
        cout << "ERROR: out of bound on real_array" << endl;;
    }
    
    tempSum[0] = 0;
    
    int error = clEnqueueWriteBuffer(commands, int_array_buff, CL_TRUE, 0, sizeof(int) *(8 + MAX_TYPE1_M_ATYPES_DIM + MAX_TYPE2_M_ATYPES_DIM + MAX_TYPE1_M_GROUPS_DIM + MAX_TYPE2_M_GROUPS_DIM) , int_array, 0, NULL, NULL);
    if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
    }
    
    error = clEnqueueWriteBuffer(commands, real_array_buff, CL_TRUE, 0, sizeof(REAL) * ( 9 + 3 + 3 + 3 + 1 + 1 + 9 + 12*3 + 3 + 9 + 12*3 + 3 + 12*3 + 12*3 + MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE + MAX_TYPE1_M_CHARGES_DIM + MAX_TYPE2_M_CHARGES_DIM), real_array, 0, NULL, NULL);
    if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
    }
    /*
     error = clEnqueueWriteBuffer(commands, temp_sum_buff, CL_TRUE, 0, sizeof(REAL), tempSum, 0, NULL, NULL);
     if(error != CL_SUCCESS){
     printf("error while writing the buffer!!\n");
     //return EXIT_FAILURE;
     }
     *//*
        error = clEnqueueWriteBuffer(commands, type1_m_atypes_dim_buff, CL_TRUE, 0, sizeof(int), &SIDECHAIN::m_aalist[getType()]->m_size, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        }
        
        error = clEnqueueWriteBuffer(commands, type2_m_atypes_dim_buff, CL_TRUE, 0, sizeof(int), &SIDECHAIN::m_aalist[pLeaf->getType()]->m_size, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        }
        
        error = clEnqueueWriteBuffer(commands, type1_m_groups_dim_buff, CL_TRUE, 0, sizeof(int), &SIDECHAIN::m_aalist[this->getType()]->m_nGroups, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        }
        
        error = clEnqueueWriteBuffer(commands, type2_m_groups_dim_buff, CL_TRUE, 0, sizeof(int), &SIDECHAIN::m_aalist[pLeaf->getType()]->m_nGroups, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        }
        
        error = clEnqueueWriteBuffer(commands, type1_m_charges_dim_buff, CL_TRUE, 0, sizeof(int), &SIDECHAIN::m_aalist[getType()]->m_size, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        }
        
        error = clEnqueueWriteBuffer(commands, type2_m_charges_dim_buff, CL_TRUE, 0, sizeof(int), &SIDECHAIN::m_aalist[pLeaf->getType()]->m_size, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        }*/
    
   	error = clEnqueueWriteBuffer(commands, type1_m_charges_buff, CL_TRUE, 0, sizeof(REAL ) * SIDECHAIN::m_aalist[getType()]->m_size,SIDECHAIN::m_aalist[this->getType()]->m_charges, 0, NULL, NULL);
        if(error != CL_SUCCESS){
        printf("error while writing the buffer!!\n");
        //return EXIT_FAILURE;
        } 
    
    //arguments for the kernel
    
    error = 0;
    error = clSetKernelArg(kernel, 0, sizeof(cl_mem), &int_array_buff);
    if (error != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments ! %d\n", error);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }
    
    error = clSetKernelArg(kernel, 1, sizeof(cl_mem), &real_array_buff);
    if (error != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments! %d\n", error);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }
    
    error = clSetKernelArg(kernel, 2, sizeof(cl_mem), &temp_sum_buff);
    if (error != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments! %d\n", error);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }
    /*
     error = clSetKernelArg(kernel, 3, sizeof(cl_mem), &type1_m_atypes_dim_buff);
     if (error != CL_SUCCESS){
     printf("Error: Failed to set kernel arguments ! %d\n", error);
     printf("Test failed\n");
     //return EXIT_FAILURE;
     }
     
     error = clSetKernelArg(kernel, 4, sizeof(cl_mem), &type2_m_atypes_dim_buff);
     if (error != CL_SUCCESS){
     printf("Error: Failed to set kernel arguments ! %d\n", error);
     printf("Test failed\n");
     //return EXIT_FAILURE;
     }
     
     error = clSetKernelArg(kernel, 5, sizeof(cl_mem), &type1_m_groups_dim_buff);
     if (error != CL_SUCCESS){
     printf("Error: Failed to set kernel arguments ! %d\n", error);
     printf("Test failed\n");
     //return EXIT_FAILURE;
     }
     
     error = clSetKernelArg(kernel, 6, sizeof(cl_mem), &type2_m_groups_dim_buff);
     if (error != CL_SUCCESS){
     printf("Error: Failed to set kernel arguments ! %d\n", error);
     printf("Test failed\n");
     //return EXIT_FAILURE;
     }
     
     error = clSetKernelArg(kernel, 7, sizeof(cl_mem), &type1_m_charges_dim_buff);
     if (error != CL_SUCCESS){
     printf("Error: Failed to set kernel arguments ! %d\n", error);
     printf("Test failed\n");
     //return EXIT_FAILURE;
     }
     
     error = clSetKernelArg(kernel, 8, sizeof(cl_mem), &type2_m_charges_dim_buff);
     if (error != CL_SUCCESS){
     printf("Error: Failed to set kernel arguments ! %d\n", error);
     printf("Test failed\n");
     //return EXIT_FAILURE;
     }
     */
    error = clSetKernelArg(kernel, 3, sizeof(cl_mem), &interaction_buff);
    if (error != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments ! %d\n", error);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }
    
    error = clSetKernelArg(kernel, 4, sizeof(cl_mem), &resetTerm_buff);
    if (error != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments ! %d\n", error);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }

    error = clSetKernelArg(kernel, 5, sizeof(cl_mem), &type1_m_charges_buff);
    if (error != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments ! %d\n", error);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }
    
    
    //computePairEnergy_k(int_array,real_array, tempSum, &SIDECHAIN::m_aalist[getType()]->m_size, &SIDECHAIN::m_aalist[pLeaf->getType()]->m_size, &SIDECHAIN::m_aalist[this->getType()]->m_nGroups, &SIDECHAIN::m_aalist[pLeaf->getType()]->m_nGroups, &SIDECHAIN::m_aalist[getType()]->m_size, &SIDECHAIN::m_aalist[pLeaf->getType()]->m_size, interaction,resetTerm);
    
    int err = 0;
    cl_event enqueue_kernel;
    err = clEnqueueTask(commands, kernel, 0, NULL, &enqueue_kernel);

    
    if (err)
    {
        printf("Error: Failed to execute kernel! %d\n", err);
        printf("Test failed\n");
        //return EXIT_FAILURE;
    }
    
    tempSum[0] = 0;
    //read results so that they are usable
    cl_event readTempSum;
    error =  clEnqueueReadBuffer(commands, temp_sum_buff, CL_TRUE, 0, sizeof(REAL), tempSum, 0, NULL, &readTempSum);
    if(error != CL_SUCCESS){
        printf("error in reading the output!! %d \n", error);
        fflush(stdout);
        //return EXIT_FAILURE;
    }
    
    cl_event readInteraction;
    error =  clEnqueueReadBuffer(commands, interaction_buff, CL_TRUE, 0, sizeof(int), interaction, 0, NULL, &readInteraction);
    if(error != CL_SUCCESS){
        printf("error in reading the output!! %d \n", error);
        fflush(stdout);
        //return EXIT_FAILURE;
    }
    
    cl_event readResetTerm;
    error =  clEnqueueReadBuffer(commands, resetTerm_buff, CL_TRUE, 0, sizeof(short), resetTerm, 0, NULL, &readResetTerm);
    if(error != CL_SUCCESS){
        printf("error in reading the output!! %d \n", error);
        fflush(stdout);
        //return EXIT_FAILURE;
    }
    
    clWaitForEvents(1, &enqueue_kernel);
    clWaitForEvents(1, &readTempSum);
    clWaitForEvents(1, &readInteraction);
    clWaitForEvents(1, &readResetTerm);
    
    
    //printf("READ: %f from sum \n", tempSum[0]);
    if (interaction[0] == 1){
        if (resetTerm[0] == 1) {
            term->reset();
        }else{
            term->set(tempSum[0]);
            m_undoPairs.push_back(term);
        }
        
    }
    
    return;
}


/*
 
 
 void CLeaf::computePairEnergy(CNode * pNode, const REAL rot[3][3],
 const REAL trans[3], CTerm * term,
 bool bSeparated)
 {
 
 assert(pNode->isLeaf());
 CLeaf * pLeaf = (CLeaf*) pNode;
 
 std::cout << pLeaf->getNext() << "\n";
 
 // GLY node does not have any atoms.
 // There is no interaction with it.
 if (getType() == GLY || pLeaf->getType() == GLY)
 return;
 
 // If the BVs are too far away, no need to do anything
 if (getBV()->computeDistance(pLeaf->getBV(), rot, trans) > CUTOFF_DISTANCE)
 {
 term->reset();
 return;
 }
 
 // Fill the pairwise distances matrix.
 computeDistances(pLeaf, rot, trans);
 
 REAL sum = 0.0;
 
 int diff = pLeaf->getIndex() - getIndex();
 
 // Compute all vdW terms.
 sum += CTerm::computeVdW(getType(), pLeaf->getType(),
 m_distances, diff);
 
 // Compute all elctrostatic terms
 sum += CTerm::computeElectrostatics(getType(), pLeaf->getType(),
 m_distances, diff);
 
 // Compute all Solvation terms.
 sum += CTerm::computeSolvation(getType(), pLeaf->getType(),
 m_distances, diff);
 
 // Store the energy sum at the corresponding leaf of the energytree.
 assert(term);
 term->set(sum);
 
 // Save a pointer to this leaf in case we need to undo the last move.
 m_undoPairs.push_back(term);
 
 return;
 }
 */

void CLeaf::computeSelfEnergy(CTerm * term)
{
    // Insert the precomputed energy of interaction between atoms inside
    // this leaf.
    term->set(m_energy);
    m_undoPairs.push_back(term);
    
    return;
}

bool CLeaf::findPairClash(CNode * pNode, const REAL rot[3][3],
                          const REAL trans[3], bool bSeparated)
{
    assert(pNode->isLeaf());
    CLeaf * pLeaf = (CLeaf*) pNode;
    
    // GLY node does not have any atoms.
    // There cannot be a clash with it.
    if (getType() == GLY || pLeaf->getType() == GLY)
        return false;
    
    int size1 = SIDECHAIN::m_aalist[getType()]->m_size;
    int size2 = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
    
    REAL cen[3], dist[3];
    
    // Check all pairs of atoms for possible clash
    for(int j = 0; j < size2; j++)
    {
        MxVpV(cen, rot, pLeaf->getPositions()[j], trans);
        int b = SIDECHAIN::m_aalist[pLeaf->getType()]->m_aTypes[j];
        
        for (int i = 0; i < size1; i++)
        {
            VmV(dist, cen, getPositions()[i]);
            int a = SIDECHAIN::m_aalist[getType()]->m_aTypes[i];
            
            // Check exclusion list to see if this pair of atoms
            // should not be checked.
            int ex = isExcluded(pLeaf->getIndex() - getIndex(), getType(), i,
                                pLeaf->getType(), j);
            if (ex != EXCLUDED &&
                isStericClash(a, b, Vlength2(dist), ex == PAIR1_4))
            {
                return true;
            }
        }
    }
    
    return false;
}

// Undo the changes to this leaf caused by the latest move.
void CLeaf::undo()
{
    if (m_nc & CNode::TRANSFORM)
    {
        m_angle = m_undoAngle;
        m_torsionE = m_undoTorsionE;
        McM(m_rotate, m_undoRotate);
    }
    
    if (m_nc & CNode::BOX)
        changeRotamer(m_undoRotIndex);
}

// Rotate the protein around the rotatable bond associated with
// this leaf.
void CLeaf::rotate(REAL angle)
{
    // Store undo information.
    m_undoAngle = m_angle;
    McM(m_undoRotate, m_rotate);
    m_undoTorsionE = m_torsionE;
    
    // Change the angle and compute a new rotation matrix.
    m_angle += angle;
    compute_rotation(m_rotate, getJoint(), m_angle);
    
    // Recompute the torsion energy for this torsion angle.
    m_torsionE = compute_dihedral(getAngle(), getIndex() % 2);
}

