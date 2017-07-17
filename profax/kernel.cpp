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


#include "kernel.h"
#include <math.h>

extern "C"{


  REAL epsilon_k[NUM_ATYPES] = {-0.1200, -0.1200, -0.0486, -0.1142, -0.1811,
			    -0.1200, -0.2384, -0.2384, -0.2384, -0.2384,
			    -0.2384, -0.2384, -0.1591, -0.1591, -0.6469,
			    -0.0430, -0.0430, -0.0498, -0.0498};

  REAL SIGMA_k[NUM_ATYPES] = {2.100, 2.100, 2.365, 2.235, 2.165, 2.100,
			  1.6000, 1.6000, 1.6000, 1.6000, 1.6000, 1.6000,
			  1.6000, 1.6000, 1.6000,  1.890, 1.890, 0.8000,
			  0.6000};

  REAL DIELECTRIC_k = 332.05382;

  REAL CUTOFF_DISTANCE_k = 9.0;

  REAL CUTOFF_DISTANCE_2_k = 9.0 * 9.0;

  REAL SOLVATION_K = 2.0 / (4.0 * M_PI * sqrt(M_PI));

  REAL deltaG_free_k[NUM_HEAVY_TYPES] = {0.00, -1.40, -0.25, 0.52, 1.50, 0.08,
				-8.90, -4.00, -7.80, -20.00, -10.00,
				-1.55, -6.70, -5.85, -10.00, -4.10,
				-2.70};

  REAL volume_k[NUM_ATYPES] = {14.7, 8.3, 23.7, 22.4, 30.0, 18.4, 4.4, 4.4,
			   11.2, 11.2, 11.2, 0.0, 10.8, 10.8, 10.8, 14.7, 21.4,
			   0.0, 0.0};

REAL
Vlength2(const REAL* V)
{
    return (V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
}

void
MxVpV(REAL *Vr,  REAL *M1,  const REAL * V1,  const REAL * V2)
{
	// REAL (* M1)[3] = ( REAL (*)[3])M1;
    	//M1 = M1;
    Vr[0] = (M1[0] * V1[0] +
             M1[1] * V1[1] +
             M1[2] * V1[2] +
             V2[0]);
    Vr[1] = (M1[3] * V1[0] +
             M1[4] * V1[1] +
             M1[5] * V1[2] +
             V2[1]);
    Vr[2] = (M1[6] * V1[0] +
             M1[7] * V1[1] +
             M1[8] * V1[2] +
             V2[2]);

}



void VmV(REAL *Vr,const REAL * V1, const REAL *  V2)
{
    Vr[0] = V1[0] - V2[0];
    Vr[1] = V1[1] - V2[1];
    Vr[2] = V1[2] - V2[2];

}

REAL Vlength(const REAL *V)
{
    return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
}

REAL computeDistance( const REAL * bv1_m_center,  const REAL * bv2_m_center,  const REAL bv1_m_rad,  const REAL bv2_m_rad,  REAL *  rot,  const REAL * trans)
{
#pragma HLS INLINE recursive

    REAL T[3];
    REAL Ttemp[3];
    
    VmV(Ttemp,trans,bv1_m_center);
    MxVpV(T, rot, bv2_m_center, Ttemp);
    
	//return 1.0;
    return Vlength(T) - (bv1_m_rad + bv2_m_rad);
}


REAL compute_vdW(int type1, int type2, REAL dist, bool is14)
{
#pragma HLS INLINE recursive
#pragma HLS PIPELINE
    // We assume "dist" is the distance squared.
    if (dist > CUTOFF_DISTANCE_2_k)
        return 0.0;
    
    REAL sig1 = SIGMA_k[type1];
    REAL eps1 = epsilon_k[type1];
    REAL sig2 = SIGMA_k[type2];
    REAL eps2 = epsilon_k[type2];
    
    if (is14)
    {
        if (type1 < NH1)
        {
            sig1 = 1.9;
            eps1 = -0.1;
        }
        
        if (type2 < NH1)
        {
            sig2 = 1.9;
            eps2 = -0.1;
        }
    }
    
    REAL eps = sqrt(eps1*eps2);
    REAL sig = (sig1 + sig2);
    
    REAL rat2 = (sig*sig)/dist;
    REAL rat6 = rat2*rat2*rat2;
    
    REAL res = (eps*(rat6*rat6 - 2*rat6));
    return res;
}

int to_index(int diff, int t1, int i, int t2, int j)
{
#pragma HLS INLINE recursive
    return (diff*1024*1024 + t1*1024*32 + i*1024 + t2*32 + j);
}

int isExcluded(int diff, int type1, int id1, int type2, int id2)
{
#pragma HLS INLINE recursive
    if (diff > 3)
        return NOT_EXCLUDED;
    else if (diff == 0 && id1 >= id2)
        return EXCLUDED;
    int index = to_index(diff,type1,id1,type2,id2);
    //case is less than 3
    if( index == 1 || index == 2 || index == 1026 || index == 1027 || index == 2051 || index == 2052 || index == 2053 || index == 3076 || index == 3077 || index == 3078 || index == 3081 || index == 5126 || index == 5129 || index == 5127 || index == 5128 || index == 5130 || index == 5131 || index == 6151 || index == 6152 || index == 6153 || index == 7176 || index == 9226 || index == 9227 || index == 10251 || index == 32801 || index == 32802 || index == 32803 || index == 33826 || index == 33827 || index == 33828 || index == 33829 || index == 35876 || index == 35877 || index == 36901 || index == 65601 || index == 65602 || index == 65603 || index == 66626 || index == 66627 || index == 67651 || index == 98401 || index == 131201 || index == 131202 || index == 132226 || index == 132227 || index == 132228 || index == 133251 || index == 133252 || index == 134276 || index == 133253 || index == 133254 || index == 134276 || index == 135301 || index == 135302 || index == 136326 || index == 164001 || index == 164002 || index == 165026 || index == 165027 || index == 165028 || index == 166051 || index == 166052 || index == 167076 || index == 196801 || index == 196802 || index == 196804 || index == 197826 || index == 197828 || index == 197827 || index == 197829 || index == 197830 || index == 198851 || index == 198854 || index == 198852 || index == 198853 || index == 199878 || index == 200901 || index == 200902 || index == 201926 || index == 229601 || index == 229602 || index == 229603 || index == 230626 || index == 231651 || index == 262401 || index == 262402 || index == 262403 || index == 263426 || index == 263427 || index == 264451 || index == 295201 || index == 295202 || index == 296226 || index == 296227 || index == 297251 || index == 297252 || index == 298276 || index == 298277 || index == 298278 || index == 298279 || index == 299301 || index == 299302 || index == 299303 || index == 300326 || index == 300327 || index == 301351 || index == 328001 || index == 328002 || index == 329026 || index == 329027 || index == 330051 || index == 360801 || index == 360802 || index == 360803 || index == 361826 || index == 361827 || index == 361828 || index == 361829 || index == 362852 || index == 362851 || index == 362854 || index == 363877 || index == 363878 || index == 364902 || index == 364901 || index == 365926 || index == 393601 || index == 393602 || index == 394626 || index == 426401 || index == 426402 || index == 427426 || index == 459201 || index == 459203 || index == 459202 || index == 460226 || index == 460227 || index == 492001 || index == 492002 || index == 492005 || index == 493026 || index == 493029 || index == 493027 || index == 493028 || index == 493030 || index == 494051 || index == 494052 || index == 494053 || index == 494054 || index == 494056 || index == 494057 || index == 495078 || index == 495080 || index == 495076 || index == 495077 || index == 495079 || index == 495082 || index == 496105 || index == 496106 || index == 497126 || index == 497127 || index == 498151 || index == 498152 || index == 500202 || index == 500201 || index == 501226 || index == 524801 || index == 524802 || index == 524804 || index == 525826 || index == 525828 || index == 525827 || index == 525829 || index == 526851 || index == 526852 || index == 526854 || index == 527878 || index == 527877 || index == 527879 || index == 528901 || index == 528902 || index == 529926 || index == 529927 || index == 530951 || index == 530952 || index == 531976 || index == 557601 || index == 557602 || index == 558626 || index == 688801 || index == 688802 || index == 689826 || index == 656001 || index == 656002 || index == 656003 || index == 656004 || index == 657026 || index == 657027 || index == 657028 || index == 658051 || index == 658052 || index == 659076 || index == 721601 || index == 721602 || index == 721603 || index == 721604 || index == 722626 || index == 723651 || index == 723652 || index == 724676 || index == 754401 || index == 754402 || index == 754403 || index == 754404 || index == 755426 || index == 756451 || index == 756452 || index == 757476 || index == 787201 || index == 787202 || index == 787203 || index == 787204 || index == 788226 || index == 789251 || index == 789252 || index == 790276 || index == 1772128 || index == 1771520 || index == 1771552 || index == 1771584 || index == 1771616 || index == 1771648 || index == 1771680 || index == 1771712 || index == 1771744 || index == 1771776 || index == 1771808 || index == 1771840 || index == 1771872 || index == 1771936 || index == 1771968 || index == 1772000 || index == 1772032 || index == 1772064 || index == 1774176 || index == 1773568 || index == 1773600 || index == 1773632 || index == 1773664 || index == 1773696 || index == 1773728 || index == 1773760 || index == 1773792 || index == 1773824 || index == 1773856 || index == 1773888 || index == 1773920 || index == 1773984 || index == 1774016 || index == 1774048 || index == 1774080 || index == 1774112 || index == 1773569 || index == 1773601 || index == 1773633 || index == 1773665 || index == 1773697 || index == 1773729 || index == 1773761 || index == 1773793 || index == 1773794 || index == 1773825 || index == 1773857 || index == 1773889 || index == 1773921 || index == 1773985 || index == 1774017 || index == 1774019 || index == 1774049 || index == 1774081 || index == 1774113 || index == 1774114 || index == 1802626 || index == 1804674 || index == 1804672 || index == 1804673 || index == 1805696 || index == 1805697 || index == 1805698 || index == 1806720 || index == 1806721 || index == 1806722 || index == 1671872 || index == 1049280 || index == 1082048 || index == 1114816 || index == 1147584 || index == 1180352 || index == 1213120 || index == 1245888 || index == 1278656 || index == 1311424 || index == 1344192 || index == 1376960 || index == 1409728 || index == 1442496 || index == 1475264 || index == 1508032 || index == 1540800 || index == 1573568 || index == 1606336 || index == 1671936 || index == 1049344 || index == 1082112 || index == 1114880 || index == 1147648 || index == 1180416 || index == 1213184 || index == 1245952 || index == 1278720 || index == 1311488 || index == 1344256 || index == 1377024 || index == 1409792 || index == 1442560 || index == 1475328 || index == 1508096 || index == 1540864 || index == 1573632 || index == 1606400 || index == 1671904 || index == 1049312 || index == 1082080 || index == 1114848 || index == 1147616 || index == 1180384 || index == 1213152 || index == 1245920 || index == 1278688 || index == 1311456 || index == 1344224 || index == 1376992 || index == 1409760 || index == 1442528 || index == 1475296 || index == 1508064 || index == 1540832 || index == 1573600 || index == 1606368 || index == 2820800 || index == 2820832 || index == 2820864 || index == 2853568 || index == 2853600 || index == 2853632 || index == 2886336 || index == 2886368 || index == 2886400 || index == 2822848 || index == 2822880 || index == 2822912 || index == 2854592 || index == 2854624 || index == 2854656 || index == 2888384 || index == 2888416 || index == 2888448 || index == 2822849 || index == 2822881 || index == 2822913 || index == 2854593 || index == 2854625 || index == 2854657 || index == 2888385 || index == 2888417 || index == 2888449 || index == 2822850 || index == 2822882 || index == 2822914 || index == 2854594 || index == 2854626 || index == 2854658 || index == 2888386 || index == 2888418 || index == 2888450 || index == 1704544 || index == 1703936 || index == 1703968 || index == 1704000 || index == 1704032 || index == 1704064 || index == 1704096 || index == 1704128 || index == 1704160 || index == 1704192 || index == 1704224 || index == 1704256 || index == 1704288 || index == 1704352 || index == 1704384 || index == 1704416 || index == 1704448 || index == 1704480 || index == 1708640 || index == 1708032 || index == 1708064 || index == 1708096 || index == 1708128 || index == 1708160 || index == 1708192 || index == 1708224 || index == 1708256 || index == 1708288 || index == 1708320 || index == 1708352 || index == 1708384 || index == 1708448 || index == 1708480 || index == 1708512 || index == 1708544 || index == 1708576 || index == 1708033 || index == 1708065 || index == 1708097 || index == 1708129 || index == 1708161 || index == 1708193 || index == 1708225 || index == 1708257 || index == 1708258 || index == 1708289 || index == 1708321 || index == 1708353 || index == 1708385 || index == 1708449 || index == 1708481 || index == 1708483 || index == 1708513 || index == 1708545 || index == 1708577 || index == 1708578 || index == 1671840 || index == 1049248 || index == 1082016 || index == 1114784 || index == 1147552 || index == 1180320 || index == 1213088 || index == 1245856 || index == 1278624 || index == 1311392 || index == 1344160 || index == 1376928 || index == 1409696 || index == 1442464 || index == 1475232 || index == 1508000 || index == 1540768 || index == 1573536 || index == 1606304 || index == 1671841 || index == 2753216 || index == 2753248 || index == 2753280 || index == 2820768 || index == 2853536 || index == 2886304 || index == 2757312 || index == 2757344 || index == 2757376 || index == 2822816 || index == 2854560 || index == 2888352 || index == 2757313 || index == 2757345 || index == 2757377 || index == 2822817 || index == 2854561 || index == 2888353 || index == 2822818 || index == 2854562 || index == 2888354 || index == 2757314 || index == 2757346 || index == 2757378){
        return EXCLUDED;
    }else if(index == 3 || index == 1028 || index == 1029 || index == 2054 || index == 2057 || index == 3079 || index == 3080 || index == 3082 || index == 3083 || index == 6154 || index == 6155 || index == 7177 || index == 8201 || index == 32804 || index == 32805 || index == 131203 || index == 131204 || index == 132229 || index == 132230 || index == 134277 || index == 134278 || index == 164003 || index == 164004 || index == 196803 || index == 196805 || index == 199876 || index == 199877 || index == 230627 || index == 295203 || index == 296228 || index == 297253 || index == 297254 || index == 297255 || index == 328003 || index == 360804 || index == 360805 || index == 361830 || index == 362853 || index == 363876 || index == 461251 || index == 492003 || index == 492004 || index == 492006 || index == 493031 || index == 493032 || index == 493033 || index == 494055 || index == 494058 || index == 495081 || index == 496101 || index == 496102 || index == 496104 || index == 497128 || index == 498154 || index == 499176 || index == 524803 || index == 524805 || index == 525830 || index == 526853 || index == 526855 || index == 527876 || index == 527880 || index == 528903 || index == 529928 || index == 722627 || index == 722628 || index == 755427 || index == 755428 || index == 788227 || index == 788228 || index == 1770080 || index == 1769472 || index == 1769504 || index == 1769536 || index == 1769568 || index == 1769600 || index == 1769632 || index == 1769664 || index == 1769696 || index == 1769728 || index == 1769760 || index == 1769792 || index == 1769824 || index == 1769888 || index == 1769920 || index == 1769952 || index == 1769984 || index == 1770016 || index == 1771521 || index == 1771553 || index == 1771585 || index == 1771617 || index == 1771649 || index == 1771681 || index == 1771713 || index == 1771745 || index == 1771746 || index == 1771777 || index == 1771809 || index == 1771841 || index == 1771873 || index == 1771937 || index == 1771969 || index == 1771971 || index == 1772001 || index == 1772033 || index == 1772065 || index == 1772066 || index == 1773152 || index == 1772544 || index == 1772576 || index == 1772608 || index == 1772640 || index == 1772672 || index == 1772704 || index == 1772736 || index == 1772768 || index == 1772800 || index == 1772832 || index == 1772864 || index == 1772896 || index == 1772960 || index == 1772992 || index == 1773024 || index == 1773056 || index == 1773088 || index == 1773570 || index == 1773602 || index == 1773603 || index == 1773634 || index == 1773635 || index == 1773698 || index == 1773730 || index == 1773762 || index == 1773764 || index == 1773795 || index == 1773826 || index == 1773827 || index == 1773858 || index == 1773890 || index == 1773922 || index == 1773923 || index == 1773986 || index == 1774018 || index == 1774050 || index == 1774053 || index == 1774082 || index == 1774084 || index == 1802624 || index == 1802625 || index == 1671873 || index == 1049281 || index == 1082049 || index == 1114817 || index == 1147585 || index == 1180353 || index == 1213121 || index == 1245889 || index == 1278657 || index == 1311425 || index == 1344193 || index == 1376961 || index == 1409729 || index == 1442497 || index == 1475265 || index == 1508033 || index == 1540801 || index == 1573569 || index == 1606337 || index == 1671874 || index == 1049282 || index == 1082050 || index == 1114818 || index == 1147586 || index == 1180354 || index == 1213122 || index == 1245890 || index == 1278658 || index == 1311426 || index == 1344194 || index == 1376962 || index == 1409730 || index == 1442498 || index == 1475266 || index == 1508034 || index == 1540802 || index == 1573570 || index == 1606338 || index == 1050304 || index == 1083072 || index == 1115840 || index == 1148608 || index == 1181376 || index == 1214144 || index == 1246912 || index == 1279680 || index == 1280704 || index == 1312448 || index == 1345216 || index == 1377984 || index == 1410752 || index == 1443520 || index == 1444544 || index == 1476288 || index == 1509056 || index == 1511104 || index == 1541824 || index == 1574592 || index == 1607360 || index == 1608384 || index == 1671937 || index == 1049345 || index == 1082113 || index == 1114881 || index == 1147649 || index == 1180417 || index == 1213185 || index == 1245953 || index == 1278721 || index == 1311489 || index == 1344257 || index == 1377025 || index == 1409793 || index == 1442561 || index == 1475329 || index == 1508097 || index == 1540865 || index == 1573633 || index == 1606401 || index == 1671938 || index == 1049346 || index == 1082114 || index == 1114882 || index == 1147650 || index == 1180418 || index == 1213186 || index == 1245954 || index == 1278722 || index == 1311490 || index == 1344258 || index == 1377026 || index == 1409794 || index == 1442562 || index == 1475330 || index == 1508098 || index == 1540866 || index == 1573634 || index == 1606402 || index == 1050368 || index == 1083136 || index == 1115904 || index == 1148672 || index == 1181440 || index == 1214208 || index == 1246976 || index == 1279744 || index == 1280768 || index == 1312512 || index == 1345280 || index == 1378048 || index == 1410816 || index == 1443584 || index == 1444608 || index == 1476352 || index == 1509120 || index == 1511168 || index == 1541888 || index == 1574656 || index == 1607424 || index == 1608448 || index == 1671905 || index == 1049313 || index == 1082081 || index == 1114849 || index == 1147617 || index == 1180385 || index == 1213153 || index == 1245921 || index == 1278689 || index == 1311457 || index == 1344225 || index == 1376993 || index == 1409761 || index == 1442529 || index == 1475297 || index == 1508065 || index == 1540833 || index == 1573601 || index == 1606369 || index == 1671906 || index == 1049314 || index == 1082082 || index == 1114850 || index == 1147618 || index == 1180386 || index == 1213154 || index == 1245922 || index == 1278690 || index == 1311458 || index == 1344226 || index == 1376994 || index == 1409762 || index == 1442530 || index == 1475298 || index == 1508066 || index == 1540834 || index == 1573602 || index == 1606370 || index == 1050336 || index == 1083104 || index == 1115872 || index == 1148640 || index == 1181408 || index == 1214176 || index == 1246944 || index == 1279712 || index == 1280736 || index == 1312480 || index == 1345248 || index == 1378016 || index == 1410784 || index == 1443552 || index == 1444576 || index == 1476320 || index == 1509088 || index == 1511136 || index == 1541856 || index == 1574624 || index == 1607392 || index == 1608416 || index == 2818752 || index == 2818784 || index == 2818816 || index == 2851520 || index == 2855616 || index == 2851552 || index == 2855648 || index == 2851584 || index == 2855680 || index == 2884288 || index == 2884320 || index == 2884352 || index == 2820801 || index == 2820833 || index == 2820865 || index == 2853569 || index == 2853601 || index == 2853633 || index == 2886337 || index == 2886369 || index == 2886401 || index == 2820802 || index == 2820834 || index == 2820866 || index == 2853570 || index == 2853602 || index == 2853634 || index == 2886338 || index == 2886370 || index == 2886402 || index == 2821824 || index == 2821856 || index == 2821888 || index == 2887360 || index == 2887392 || index == 2887424 || index == 2822851 || index == 2822915 || index == 2854595 || index == 2854659 || index == 2888387 || index == 2888451 || index == 2822852 || index == 2822883 || index == 2822884 || index == 2822916 || index == 2854596 || index == 2854627 || index == 2854628 || index == 2854660 || index == 2888388 || index == 2888419 || index == 2888419 || index == 2888452 || index == 1703937 || index == 1703969 || index == 1704001 || index == 1704033 || index == 1704065 || index == 1704097 || index == 1704129 || index == 1704161 || index == 1704162 || index == 1704193 || index == 1704225 || index == 1704257 || index == 1704289 || index == 1704353 || index == 1704385 || index == 1704387 || index == 1704417 || index == 1704449 || index == 1704481 || index == 1704482 || index == 1705568 || index == 1704960 || index == 1704992 || index == 1705024 || index == 1705056 || index == 1705088 || index == 1705120 || index == 1705152 || index == 1705184 || index == 1705216 || index == 1705248 || index == 1705280 || index == 1705312 || index == 1705376 || index == 1705408 || index == 1705440 || index == 1705472 || index == 1705504 || index == 1706592 || index == 1705984 || index == 1706016 || index == 1706048 || index == 1706080 || index == 1706112 || index == 1706144 || index == 1706176 || index == 1706208 || index == 1706240 || index == 1706272 || index == 1706304 || index == 1706336 || index == 1706400 || index == 1706432 || index == 1706464 || index == 1706496 || index == 1706528 || index == 1707616 || index == 1707008 || index == 1707040 || index == 1707072 || index == 1707104 || index == 1707136 || index == 1707168 || index == 1707200 || index == 1707232 || index == 1707264 || index == 1707296 || index == 1707328 || index == 1707360 || index == 1707424 || index == 1707456 || index == 1707488 || index == 1707520 || index == 1707552 || index == 1708034 || index == 1708066 || index == 1708067 || index == 1708098 || index == 1708099 || index == 1708162 || index == 1708194 || index == 1708226 || index == 1708228 || index == 1708259 || index == 1708290 || index == 1708291 || index == 1708322 || index == 1708354 || index == 1708386 || index == 1708387 || index == 1708450 || index == 1708482 || index == 1708514 || index == 1708517 || index == 1708546 || index == 1708548 || index == 1049249 || index == 1082017 || index == 1114785 || index == 1147553 || index == 1180321 || index == 1213089 || index == 1245857 || index == 1278625 || index == 1311393 || index == 1344161 || index == 1376929 || index == 1409697 || index == 1442465 || index == 1475233 || index == 1508001 || index == 1540769 || index == 1573537 || index == 1606305 || index == 1671842 || index == 1049250 || index == 1082018 || index == 1114786 || index == 1147554 || index == 1180322 || index == 1213090 || index == 1245858 || index == 1278626 || index == 1311394 || index == 1344162 || index == 1376930 || index == 1409698 || index == 1442466 || index == 1475234 || index == 1508002 || index == 1540770 || index == 1573538 || index == 1606306 || index == 1050272 || index == 1083040 || index == 1115808 || index == 1148576 || index == 1181344 || index == 1214112 || index == 1246880 || index == 1279648 || index == 1312416 || index == 1345184 || index == 1377952 || index == 1410720 || index == 1443488 || index == 1476256 || index == 1509024 || index == 1541792 || index == 1574560 || index == 1607328 || index == 2753217 || index == 2753249 || index == 2753281 || index == 2820769 || index == 2853537 || index == 2886305 || index == 2753218 || index == 2753250 || index == 2753282 || index == 2820770 || index == 2853538 || index == 2886306 || index == 2754240 || index == 2754272 || index == 2754304 || index == 2755264 || index == 2755296 || index == 2755328 || index == 2756288 || index == 2756320 || index == 2756352 || index == 2818720 || index == 2851488 || index == 2855584 || index == 2884256 || index == 2757315 || index == 2757379 || index == 2757316 || index == 2757347 || index == 2757348 || index == 2757380 || index == 3871106 || index == 3936642 || index == 3902850 || index == 3805570){
        return PAIR1_4;
    }else{
        return NOT_EXCLUDED;
    }
    /*
     EXCLUSIONS::const_iterator ex = exclusion_list.find(to_index(diff,type1,id1,type2,id2));
     
     if (ex == exclusion_list.end())
     return NOT_EXCLUDED;
     else if (ex->second < 3)
     return EXCLUDED;
     else if (ex->second == 3)
     return PAIR1_4;
     else
     assert(false);
     */
    //return 0;
}

REAL computeVdW(int size1, int size2, int type1, int type2,  int * type1_m_aTypes,  int * type2_m_aTypes, REAL * dists,int diff)
{
#pragma HLS INLINE recursive
    //assert(diff >= 0);
    
    REAL sum[1];
    sum[0]= 0.0;
    
    compute_vdw_outer:for (int i = 0; i < size1; i++)
    {

    	compute_vdw_inner:for (int j = 0; j < size2; j++)
        {
#pragma HLS PIPELINE
            int ex = isExcluded(diff, type1, i, type2, j);
            if (ex != EXCLUDED)
            {
                sum[0] += compute_vdW(type1_m_aTypes[i],
                                   type2_m_aTypes[j],
                                   dists[i*MAX_ROTAMER_SIZE + j], ex == PAIR1_4);
            }
        }
    }
    return sum[0];
}

REAL compute_ES(REAL charge1, REAL charge2, REAL dist)
{
#pragma HLS INLINE recursive
    // We assume the "dist" is the squared distance.
    return (DIELECTRIC_k * (charge1 * charge2)/(dist));
}


REAL computeElectrostatics(int size1, int size2, int type1_m_nGroups, int type2_m_nGroups,  int * type1_m_groups, int * type2_m_groups, REAL * type1_m_charges,  REAL * type2_m_charges, int type1, int type2,
                            REAL * dists,
                           int diff)
{
#pragma HLS INLINE recursive
    //assert(diff >= 0);
    
    // Proline is treated differently because the Cd atom is in a group
    // with the Ca atom which is in a different leaf node.
    if (type1 == BBP)
        size1++;
    else if (type1 == PRO)
        size1--;
    
    if (type2 == BBP)
        size2++;
    else if (type2 == PRO)
        size2--;
    
    int ngr1 = type1_m_nGroups;
    int ngr2 = type2_m_nGroups;
    
    REAL sum = 0.0;
    
    for (int g1 = 0; g1 < ngr1; g1++)
        for (int g2 = 0; g2 < ngr2; g2++)
        {
            int end1, end2, start1, start2;
            
            start1 = type1_m_groups[g1];
            start2 = type2_m_groups[g2];
            
            if (g1 == ngr1 - 1)
                end1 = size1;
            else
                end1 = type1_m_groups[g1 + 1];
            
            if (g2 == ngr2 - 1)
                end2 = size2;
            else
                end2 = type2_m_groups[g2 + 1];
            
            // Decide whether the two groups are close enough to interact.
            // They are close enough if they contain an interacting pair
            // of atoms.

            bool bCompute = false;
            ce_loop1:for (int i = start1; i < end1; i++)
                c2_loop2:for(int j = start2; j < end2 && i!=1000; j++)
                {
        #pragma HLS PIPELINE
                    if (dists[i*MAX_ROTAMER_SIZE + j] < CUTOFF_DISTANCE_2_k)
                    {
                        bCompute = true;
                        i = 1000;
                        //break;
                    }
                }
            
            // If the groups are interacting, compute their contribution.
            if (bCompute)
            {
                for (int i = start1; i < end1; i++)
                    for(int j = start2; j < end2; j++)
                    {
                        int ex = isExcluded(diff, type1, i, type2, j);
                        if (ex != EXCLUDED)
                        {
                            REAL w = compute_ES(type1_m_charges[i], type2_m_charges[j],
                                                dists[i*MAX_ROTAMER_SIZE + j]);
                            if (ex == PAIR1_4)
                                w *= 0.4;
                            sum += w;
                        }
                    }
            }
        }
    
    return sum;
}

REAL computeSolventEffect(int aType1, int aType2, REAL dist, REAL lambda1, REAL lambda2, bool is14)
{
#pragma HLS INLINE recursive
    if (dist > CUTOFF_DISTANCE_2_k)
        return 0.0;
    
    REAL t1 = 0.0, t2 = 0.0;
    
    if (aType1 < NUM_HEAVY_TYPES && aType2 < NUM_HEAVY_TYPES)
    {
        REAL sig1 = SIGMA_k[aType1], sig2 = SIGMA_k[aType2];
        
        if (is14)
        {
            if (aType1 < NH1)
                sig1 = 1.9;
            
            if (aType2 < NH1)
                sig2 = 1.9;
        }
        
        // We assume "dist" is the squared distance.
        REAL d = sqrt(dist);
        
        REAL X12 = (d - sig1)/lambda1;
        t1 = SOLVATION_K * deltaG_free_k[aType1] * exp(-X12*X12) * volume_k[aType2]
        / (lambda1 * dist);
        
        REAL X21 = (d - sig2)/lambda2;
        t2 = SOLVATION_K * deltaG_free_k[aType2] * exp(-X21*X21) * volume_k[aType1]
        / (lambda2 * dist);
    }
    
    return -(t2 + t1);
}

REAL getLambda(int AAtype, int index, int aType)
{
#pragma HLS INLINE recursive
    if (aType == NH3 || aType == NC2 || aType == OC)
        return 6.0;
    
    if ((AAtype == ARG && index >= 2) ||
        (AAtype == LYS && index >= 3) ||
        (AAtype == ASP) ||
        (AAtype == GLU && index >= 1) ||
        (AAtype == NTR) ||
        (AAtype = CTR))
        return 6.00;
    else
        return 3.50;
}

REAL computeSolvation(int size1, int size2,  int * type1_m_aTypes, int * type2_m_aTypes, int type1, int type2,
                       REAL *  dists,
                      int diff)
{
#pragma HLS INLINE recursive
    //assert(diff >= 0);
    
    REAL sum = 0.0;
    REAL tempSum[16]={0};
    
    REAL tempSum8[8];
  REAL tempSum4[4];
  REAL tempSum2[2];

solvation_out:for (int i = 0; i < 12; i++) {
  #pragma HLS PIPELINE II=1

  for(int i=0; i < 16; i++) {
    tempSum[i]=0;
  }
      if(i<size1){
        for (int j = 0; j < 12; j++)
        {
          if(j<size2) {
            int ex = isExcluded(diff, type1, i, type2, j);
            if (ex != EXCLUDED)
            {
                REAL lambda1 = getLambda(type1, i, type1_m_aTypes[i]);
                REAL lambda2 = getLambda(type2, j, type2_m_aTypes[j]);
    REAL dists_param = dists[i*MAX_ROTAMER_SIZE + j];
        /*
        sum += computeSolventEffect(type1_m_aTypes[i],
                                            type2_m_aTypes[j],
                                            dists_param, lambda1, lambda2, ex == PAIR1_4);
        */
        tempSum[j]= computeSolventEffect(type1_m_aTypes[i],
                      type2_m_aTypes[j],
                      dists_param, lambda1, lambda2, ex == PAIR1_4);


            }
        }
        }

        for(int s=0; s<8; s++)
      tempSum8[s] = tempSum[2*s+0] + tempSum[2*s+1];

    for(int s=0; s<4; s++)
      tempSum4[s] = tempSum8[2*s+0] + tempSum8[2*s+1];

    for(int s=0; s<2; s++)
      tempSum2[s] = tempSum4[2*s+0] + tempSum4[2*s+1];

    sum += tempSum2[0] + tempSum2[1];
      }
    }
    return sum;
}


 void computeDistances(int size1, int size2, int CLeafType, int pNodeType,  REAL *  CLeaf_next_m_positions,  REAL * CLeaf_m_translate,  REAL * CLeaf_m_rotate, REAL *  pNode_next_m_positions,  REAL * pNode_m_translate,  REAL * pNode_m_rotate, REAL *  pNode_m_positions,  REAL * CLeaf_m_positions,  REAL * CLeaf_m_distances, REAL * rot,  const REAL * trans){
#pragma HLS INLINE recursive
    
    REAL cen[3], dist[3], vec1[3], vec2[3];
  
  
    // If this is a PRO backbone, we need to add the Cd atom because
    // it is part of a group with the Ca atom for Electrostatic purposes
    if (CLeafType == BBP)
    {
        //assert(getNext());
        //assert(getNext()->getType() == PRO);
        MxVpV(vec1, CLeaf_m_rotate, CLeaf_next_m_positions + 6, CLeaf_m_translate); // 6 = i = 2
    }
    
    // If this is a PRO backbone, we need to add the Cd atom because
    // it is part of a group with the Ca atom for Electrostatic purposes
    if (pNodeType == BBP)
    {
        //assert(pLeaf->getNext());
        //assert(pLeaf->getNext()->getType() == PRO);
        
        REAL temp[3];
  MxVpV(temp, pNode_m_rotate, pNode_next_m_positions + 6,
              pNode_m_translate);
        MxVpV(vec2, rot,temp, trans);
    }
    
    // Compute te distances between all pairs of atoms.

    for(int j = 0; j < 12; j++)

    {

      if(j<size2) {

  //REAL *pNode_m_positions_p[3];
  //for(int i = 0; i < 3; i++) pNode_m_positions_p[i] = TODO CHECK OTHER DIMENSION....
        //MxVpV(cen, rot, pNode_m_positions[j], trans);
        MxVpV(cen, rot, pNode_m_positions+ j * ROW, trans);
        
  for (int i = 0; i < 12; i++){
#pragma HLS PIPELINE II=1
    {
    if(i<size1) {
            VmV(dist, cen, CLeaf_m_positions + i*ROW);
        
            CLeaf_m_distances[i * MAX_ROTAMER_SIZE + j] = Vlength2(dist);
            // Add distances to the Cd of the second node (if type is BBP)
            if (pNodeType == BBP)
            {
                VmV(dist, vec2, CLeaf_m_positions + i * ROW);
                CLeaf_m_distances[i* MAX_ROTAMER_SIZE + size2] = Vlength2(dist);
            }
    }
  }
  
        
  }
        // Add distances to the Cd of the first node (if type is BBP)
        if (CLeafType == BBP)
        {
            VmV(dist, cen, vec1);
            CLeaf_m_distances[size1*MAX_ROTAMER_SIZE + j] = Vlength2(dist);
        }

     
    // Add distance between the Cd of the first and second nodes (both BBPs)
    if (CLeafType == BBP && pNodeType == BBP)
    {
        VmV(dist, vec2, vec1);
        CLeaf_m_distances[size1*MAX_ROTAMER_SIZE + size2] = Vlength2(dist);
    }
    }
    }

  
  return;    
}


void computePairEnergy_k( int * int_array,  REAL * real_array,  REAL * termSum,  int* interaction,  short *resetTerm, REAL * type1_m_charges_g){
 
#pragma HLS INTERFACE m_axi port=int_array offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=real_array offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=termSum offset=slave bundle=gmem2
#pragma HLS INTERFACE m_axi port=interaction offset=slave bundle=gmem9
#pragma HLS INTERFACE m_axi port=resetTerm offset=slave bundle=gmem10
#pragma HLS INTERFACE m_axi port=type1_m_charges_g offset=slave bundle=gmem3

#pragma HLS INTERFACE s_axilite port=int_array bundle=control
#pragma HLS INTERFACE s_axilite port=real_array bundle=control
#pragma HLS INTERFACE s_axilite port=termSum bundle=control
#pragma HLS INTERFACE s_axilite port=interaction bundle=control
#pragma HLS INTERFACE s_axilite port=resetTerm bundle=control
#pragma HLS INTERFACE s_axilite port=type1_m_charges_g bundle=control

#pragma HLS INTERFACE s_axilite port=return bundle=control

	//read the input buffer...
	int CLeafType = int_array[0];
	int pNodeType = int_array[1];
	int size1 = int_array[2];
	int size2 = int_array[3];
    int CLeaf_m_index = int_array[4];
    int pNode_m_index = int_array[5];
    
    int type1_m_aTypes[MAX_TYPE1_M_ATYPES_DIM];
    int type2_m_aTypes[MAX_TYPE2_M_ATYPES_DIM];
    
    int type1_m_groups[MAX_TYPE1_M_GROUPS_DIM];
    int type2_m_groups[MAX_TYPE2_M_GROUPS_DIM];
    
	int index = 6;

	read_t1mat:for(int i = index, s_i = 0; i < index + size1; i++, s_i++){
#pragma HLS PIPELINE II=1
    	type1_m_aTypes[s_i] = int_array[i];
    }
	//index += TYPE1_M_ATYPES_DIM;
	index += size1;
	
	read_t2mat:for(int i = index, s_i = 0; i < index + size2; i++, s_i++) {
#pragma HLS PIPELINE II=1
		type2_m_aTypes[s_i] = int_array[i];
	}
	//index += TYPE2_M_ATYPES_DIM;
	index += size2;

	
	int type1_m_nGroups = int_array[index];
	index++;
	int type2_m_nGroups = int_array[index];
	index++;

	read_t1mg:for(int i = index, s_i = 0; i < index + type1_m_nGroups; i++, s_i++){
#pragma HLS PIPELINE II=1
		type1_m_groups[s_i] = int_array[i];
	}
	//index += TYPE1_M_GROUPS_DIM;
    index += type1_m_nGroups;
	read_t2mg:for(int i = index, s_i = 0; i < index + type2_m_nGroups; i++, s_i++) {
#pragma HLS PIPELINE II=1
		type2_m_groups[s_i] = int_array[i];
	}	

	REAL rot[3 * 3];
	REAL trans[3];
	REAL bv1_m_center[3];
	REAL bv2_m_center[3];
	REAL bv1_m_rad[1];
	REAL bv2_m_rad[1];
	REAL CLeaf_m_rotate[3 * 3];
	REAL CLeaf_next_m_positions[12 * 3];
	REAL CLeaf_m_translate[3];
	REAL pNode_m_rotate[3 * 3];
	REAL pNode_next_m_positions[12 * 3];
	REAL pNode_m_translate[3];
	REAL pNode_m_positions[12 * 3];
	REAL CLeaf_m_positions[12 * 3];
	REAL CLeaf_m_distances[MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE];
	REAL type1_m_charges[MAX_TYPE1_M_CHARGES_DIM];
	REAL type2_m_charges[MAX_TYPE2_M_CHARGES_DIM];

	index = 0;
	read_rot:for(int i = index, s_i = 0; i < index + 9; i++, s_i++){
#pragma HLS PIPELINE
		rot[s_i] = real_array[i];
	}
	index += 9;
	read_trans:for(int i = index, s_i = 0; i < index + 3; i++, s_i++){
#pragma HLS PIPELINE
            trans[s_i] = real_array[i];
        }
        index += 3;
        read_bv1mc:for(int i = index, s_i = 0; i < index + 3; i++, s_i++){
#pragma HLS PIPELINE
            bv1_m_center[s_i] = real_array[i];
        }
        index +=3;
        read_bv2mc:for(int i = index, s_i = 0; i < index + 3; i++, s_i++){
#pragma HLS PIPELINE
            bv2_m_center[s_i] = real_array[i];
        }
        index +=3;
        bv1_m_rad[0] = real_array[index];
        index ++;
        bv2_m_rad[0] = real_array[index];
        index ++;
        read_cmr:for(int i = index, s_i = 0; i < index + 9; i++, s_i++){
#pragma HLS PIPELINE
            CLeaf_m_rotate[s_i] = real_array[i];
        }
        index += 9;
        read_cnmp:for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++){
#pragma HLS PIPELINE
            CLeaf_next_m_positions[s_i] = real_array[i];
        }
        index += 12 * 3;
        read_cmt:for(int i = index, s_i = 0; i < index + 3; i++, s_i++){
#pragma HLS PIPELINE
            CLeaf_m_translate [s_i] = real_array[i];
        }
        index += 3;
        read_pmr:for(int i = index, s_i = 0; i < index + 9; i++, s_i++){
#pragma HLS PIPELINE
            pNode_m_rotate[s_i] = real_array[i];
        }
        index += 9;
        read_pnmp:for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++){
#pragma HLS PIPELINE
            pNode_next_m_positions[s_i] = real_array[i];
        }
        index += 12 * 3;
        read_pmt:for(int i = index, s_i = 0; i < index + 3; i++, s_i++){
#pragma HLS PIPELINE
            pNode_m_translate[s_i] = real_array[i];
        }
        index += 3;
        read_pmnp:for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++){
#pragma HLS PIPELINE
            pNode_m_positions[s_i] = real_array[i];
        }
        index += 12 * 3;
        read_cmp:for(int i = index, s_i = 0; i < index + 12 * 3; i++, s_i++){
#pragma HLS PIPELINE
            CLeaf_m_positions[s_i] = real_array[i];
        }
        index += 12 * 3;
        read_cmd:for(int i = index, s_i = 0; i < index + MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE; i++, s_i++){
#pragma HLS PIPELINE
            CLeaf_m_distances[s_i] = real_array[i];
        }
        index += MAX_ROTAMER_SIZE * MAX_ROTAMER_SIZE;
        /*read_t1mc:for(int i = index, s_i = 0; i < index + size1; i++, s_i++){
#pragma HLS PIPELINE
            type1_m_charges[s_i] = real_array[i];
        }*/

        read_t1mc:for(int i = 0; i < size1; i++){
#pragma HLS PIPELINE
            type1_m_charges[i] = type1_m_charges_g[i];
        }
        index += size1;
        read_t2mc:for(int i = index, s_i = 0; i < index + size2; i++, s_i++){
#pragma HLS PIPELINE
            type2_m_charges[s_i] = real_array[i];
        }

        //END READING

	//for(int i = 0; i < size1; i++) type1_m_charges[i] += 0;

    if(CLeafType == GLY || pNodeType == GLY){
        interaction[0] = 0; //false
        //return; ZIO POVERO!
        *termSum = 0; //tanto non se li caga se interaction[0] è falso
        resetTerm[0] = 0; //tanto non se li caga se interaction[0] è falso
    }else{
        interaction[0] = 1; //means true
    // If the BVs are too far away, no need to do anything


        if ( computeDistance(bv1_m_center, bv2_m_center, *(bv1_m_rad),*(bv2_m_rad), rot, trans) > CUTOFF_DISTANCE_k)
        {
                resetTerm[0] = 1; //true
                //return;
                *termSum = 0;
        }else{
            resetTerm[0] = 0; //false
            computeDistances(size1, size2, CLeafType, pNodeType, CLeaf_next_m_positions, CLeaf_m_translate, CLeaf_m_rotate, pNode_next_m_positions, pNode_m_translate, pNode_m_rotate, pNode_m_positions, CLeaf_m_positions, CLeaf_m_distances, rot,  trans);


        
     
            REAL sum = 0.0;
            
            int diff = pNode_m_index - CLeaf_m_index;

            // Compute all vdW terms.
            sum += computeVdW(size1, size2, CLeafType, pNodeType, type1_m_aTypes, type2_m_aTypes,
                              CLeaf_m_distances, diff);
            // Compute all elctrostatic terms
            sum += computeElectrostatics(size1, size2, type1_m_nGroups, type2_m_nGroups, type1_m_groups, type2_m_groups, type1_m_charges, type2_m_charges, CLeafType, pNodeType, CLeaf_m_distances, diff);
            
            // Compute all Solvation terms.
            sum += computeSolvation(size1, size2, type1_m_aTypes, type2_m_aTypes, CLeafType, pNodeType,
                                    CLeaf_m_distances, diff);
            

            *termSum = sum;

        }
    

    
    }
    return;
    
    	
}

}
