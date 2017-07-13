/*
 * Copyright 2016 Lorenzo Di Tucci, Giulia Guidi, Emanuele Del Sozzo
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * 	http://www.apache.org/licenses/LICENSE-2.0
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



#ifndef _LEAFCL_HPP_
#define _LEAFCL_HPP_

#include "eef1.hpp"


void computePairEnergyCL(int CLeafType, int pNodeType, bool pNodeIsLeaf, const REAL rot[3][3],
                         const REAL trans[3], bool bSeparated, const REAL bv1_m_center[3], const REAL bv2_m_center[3], REAL bv1_m_rad, REAL bv2_m_rad, int size1, int size2, REAL CLeaf_m_rotate[3][3], COORDS CLeaf_next_m_positions, REAL CLeaf_m_translate[3], REAL pNode_m_rotate[3][3], COORDS pNode_next_m_positions, REAL pNode_m_translate[3], COORDS pNode_m_positions, COORDS CLeaf_m_positions, REAL CLeaf_m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE], int CLeaf_m_index, int pNode_m_index, int * type1_m_aTypes, int * type2_m_aTypes, int type1_m_nGroups, int type2_m_nGroups, int * type1_m_groups, int * type2_m_groups, REAL * type1_m_charges, REAL * type2_m_charges);


REAL computeDistance(const REAL bv1_m_center[3], const REAL bv2_m_center[3], REAL bv1_m_rad, REAL bv2_m_rad, const REAL rot[3][3], const REAL trans[3]);

void computeDistances(int size1, int size2, int CLeafType, int pNodeType, COORDS CLeaf_next_m_positions, REAL CLeaf_m_translate[3], REAL CLeaf_m_rotate[3][3], COORDS pNode_next_m_positions, REAL pNode_m_translate[3], REAL pNode_m_rotate[3][3], COORDS pNode_m_positions, COORDS CLeaf_m_positions, REAL CLeaf_m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE], const REAL rot[3][3], const REAL trans[3]);

REAL computeVdW(int size1, int size2, int type1, int type2, int * type1_m_aTypes, int * type2_m_aTypes,
                REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE], int diff);

REAL computeElectrostatics(int size1, int size2, int type1_m_nGroups, int type2_m_nGroups, int * type1_m_groups, int * type2_m_groups, REAL * type1_m_charges, REAL * type2_m_charges, int type1, int type2,
                           REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE],
                           int diff);

REAL computeSolvation(int size1, int size2, int * type1_m_aTypes, int * type2_m_aTypes, int type1, int type2,
                      REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE],
                      int diff);

#endif
