#include <string.h>
#include <ap_int.h>

#define NUM_ELEM_BITS 256
#define NUM_ELEM 128
#define SIZE_BYTE 32
#define N 256
//#define M 1048576
#define M 512

const short GAP_i = -1;
const short GAP_d = -1;
const short MATCH = 2;
const short MISS_MATCH = -1;
const short CENTER = 0;
const short NORTH = 1;
const short NORTH_WEST = 2;
const short WEST = 3;
extern "C" {




void update_database(ap_uint<NUM_ELEM_BITS> *database, ap_uint<NUM_ELEM_BITS> *databaseLocal, int num_diagonals){
	int startingIndex = N + num_diagonals;
	update_database:for(int i = 1; i < N; i++){
#pragma HLS PIPELINE
		databaseLocal[(i-1)/NUM_ELEM].range(((i-1)%NUM_ELEM)*2+1, ((i-1)%NUM_ELEM)*2) = databaseLocal[i/NUM_ELEM].range((i%NUM_ELEM)*2+1, (i%NUM_ELEM)*2);
		//set_char(databaseLocal, i-1, get_char(databaseLocal, i));
		//databaseLocal[i-1] = databaseLocal[i];
	}
	databaseLocal[(N-1)/NUM_ELEM].range(((N-1)%NUM_ELEM) * 2 + 1,((N-1)%NUM_ELEM) * 2) = database[startingIndex/NUM_ELEM].range((startingIndex%NUM_ELEM) * 2 +1, (startingIndex%NUM_ELEM) *2);
	//set_2bit(databaseLocal, N-1, get_2bit(database, startingIndex));
	//set_char(databaseLocal, N-1, get_char(database, startingIndex));
	//databaseLocal[N-1] = database[startingIndex];
}
  void calculate_and_store(short *north, short *northwest, short *west, ap_uint<NUM_ELEM_BITS> *query, ap_uint<NUM_ELEM_BITS> *database, int Nl, int Ml, int num_diagonals, ap_int<512> *directionMatrix, int queryIndex, int databaseIndex,ap_int<512> *temp_p, int to[1], int from[1], int *elem_stored, int *store_offset, ap_uint<NUM_ELEM_BITS> *databaseLocal){
  #pragma HLS INLINE region recursive
    //calculate_diagonal(north, northwest, west, query, database,directionDiagonal, Nl,Ml, num_diagonals, queryIndex, databaseIndex);
	  //calculate one diagonal
	  	  	  //from[0] = 510;
	  	  	  //to[0] = 511;
	  	  	  from[0] = N * 2 - 2;
	  	  	  to[0] = N * 2 - 1;
	  	  	  int databaseLocalIndex = 0;
	      	        calculate_diagonal_for:for(int index = N-1; index >=0 ; index--){

	      	  #pragma HLS PIPELINE
	      	            int val = 0;
	      	            //unsigned int q = get_2bit(query, index);
	      	            //unsigned int db = get_2bit(databaseLocal, databaseLocalIndex);
	      	            unsigned int q = query[index/NUM_ELEM].range((index%NUM_ELEM) * 2+ 1, (index % NUM_ELEM) *2);
	      	            unsigned int db = databaseLocal[databaseLocalIndex/NUM_ELEM].range((databaseLocalIndex % NUM_ELEM) * 2 + 1, (databaseLocalIndex % NUM_ELEM) * 2);
	      	            if(num_diagonals < N - 1 && databaseLocalIndex < N - 1 - num_diagonals) db = 9;
	      	            //printf("q is %d, db is %d \n", q, db);
	      	            //printf(" %d", db);
	      	            const short match = (q == db) ? MATCH:MISS_MATCH;
	      	            //const short match = (query[index] == databaseLocal[databaseLocalIndex]) ? MATCH:MISS_MATCH;
	      	            const short val1 = northwest[index] + match;
	      	            const short val2 = north[index] + GAP_d;
	      	            const short val3 = west[index] + GAP_i;

	      	            if(val1 > val && val1 >= val2 && val1 >= val3){
	      	              //val1
	      	              *(northwest + index + 1) = *(north + index);
	      	              *(north + index) = val1;
	      	              *(west + index + 1) = val1;
	      	              //*(directionDiagonal + index) = NORTH_WEST;
	      	              temp_p[0].range(to[0], from[0]) = NORTH_WEST;
	      	            }else if(val2 > val && val2 >= val3){
	      	              //val2
	      	              *(northwest + index + 1) = *(north + index);
	      	              *(north + index) = val2;
	      	              *(west + index + 1) = val2;
	      	              //*(directionDiagonal + index) = NORTH;
	      	              temp_p[0].range(to[0], from[0]) = NORTH;
	      	            }else if(val3 > val){
	      	              //val3
	      	              *(northwest + index + 1) = *(north + index);
	      	              *(north + index) = val3;
	      	              *(west + index + 1) = val3;
	      	              //*(directionDiagonal + index) = WEST;
	      	              temp_p[0].range(to[0], from[0]) = WEST;
	      	            }else{
	      	              //val
	      	              *(northwest + index + 1) = *(north + index);
	      	              *(north + index) = val;
	      	              *(west + index + 1) = val;
	      	              //*(directionDiagonal + index) = CENTER;
	      	              temp_p[0].range(to[0], from[0]) = CENTER;
	      	            }

	      	            to[0] -=2;
	      	            from[0] -=2;
	      	            //databaseIndex++;
	      	            databaseLocalIndex++;



	      	        }//endloop for calculating diagonal
	      	        //printf("\n");
	      	      	memcpy((ap_int<512>*) ((directionMatrix + num_diagonals)), temp_p,  64);


  }

  void smithwaterman(ap_uint<NUM_ELEM_BITS>  *g_query, ap_uint<NUM_ELEM_BITS> *g_database, ap_int<512> *directionMatrix)
  {
  #pragma HLS INTERFACE m_axi port=g_query offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=g_database offset=slave bundle=gmem1
  #pragma HLS INTERFACE m_axi port=directionMatrix offset=slave bundle=gmem2

  #pragma HLS INTERFACE s_axilite port=g_query bundle=control
  #pragma HLS INTERFACE s_axilite port=g_database bundle=control
  #pragma HLS INTERFACE s_axilite port=directionMatrix bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

	  ap_uint<NUM_ELEM_BITS> query[N/NUM_ELEM];
//#pragma HLS ARRAY_PARTITION variable=query complete dim=1
	  ap_uint<NUM_ELEM_BITS> database[(M + 2*(N))/NUM_ELEM];
//#pragma HLS ARRAY_PARTITION variable=database complete dim=1

	  ap_uint<NUM_ELEM_BITS> databaseLocal[N/NUM_ELEM];
#pragma HLS ARRAY_PARTITION variable=databaseLocal complete dim=1

    memcpy(query, g_query, N/NUM_ELEM * SIZE_BYTE);
    memcpy(database, g_database, (M+2*(N))/NUM_ELEM * SIZE_BYTE);

     int num_diagonals;
    
    //buffer for the directiondiagonal
  // short directionDiagonal[N];
  //#pragma HLS ARRAY_PARTITION variable=directionDiagonal complete dim=1
    //int store_offset[1];
    //store_offset[0] = 0; //offset for the directionMatrix
    //buffer needed for dependencies
    short north[N + 1];
  #pragma HLS ARRAY_PARTITION variable=north complete dim=1
    short west[N + 1];
  #pragma HLS ARRAY_PARTITION variable=west complete dim=1
    short northwest[N + 1];
  #pragma HLS ARRAY_PARTITION variable=northwest complete dim=1

    ap_int<512> temp_p[1];
    #pragma HLS ARRAY_PARTITION variable=temp_p complete dim=1
    int to[1];
    int from[1];
    int elem_stored[1];
    int store_offset[1];

    store_offset[0] = 0;

    to[0] = 1;
    from[0] = 0;
    elem_stored[0] = 0;



    //init dependecy buffers
    init_dep_for:for(int i = 0; i<= N; i++){
      #pragma HLS PIPELINE
      north[i] = 0;
      west[i] = 0;
      northwest[i] = 0;
    }
    
    //init the database with the first N values.
    init_db:for(int i = 0; i < N/NUM_ELEM; i++){
#pragma HLS PIPELINE
    	databaseLocal[i] = database[i];
    }

    //initi query and DB
    //int databaseIndex = 0;
    int queryIndex = N - 1;
    int index = 0;
      outer_for:for(num_diagonals = 0; num_diagonals < N+M-1; num_diagonals++){
   //#pragma HLS DEPENDENCE variable=temp_p false
  #pragma HLS inline region recursive
  #pragma HLS PIPELINE


            calculate_and_store(north, northwest, west, query, database, N,M, num_diagonals, directionMatrix, queryIndex, num_diagonals,temp_p, to, from, elem_stored, store_offset, databaseLocal);
            update_database(database,databaseLocal,num_diagonals);

      }//end outer for
      //store_memory_all(temp_p, directionMatrix);
      //memcpy((ap_int<512>*) directionMatrix, temp_p,  NUM_ELEM_BITS * ((N * (N+M-1) / 256) + 1));
      //*maxIndex = localMaxIndex;
      //*maxIndexX = localMaxIndexX;
      //*maxIndexY = localMaxIndexY;
  return;
  }
}
