#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "helper.h"
#define SERIAL_ELEMENT_LIMIT 11
#define NUM_BLOCKS 32
#define NUM_CHUNK 32
#define NUM_DIVIDE 64

int main(int argc, char*argv[])
{

  Matrix *A, *B, *C, *F, *D, *result;
/*
  int cola[4] = {0,2,3,5};
  int rowa[5] = {1,2,2,0,2};
  A = value_creator(3,3,rowa,cola);
  int colb[4] = {0,0,2,4};
  int rowb[4] = {0,2,1,2};
  B = value_creator(3,3,rowb,colb);
  int colc[4] = {0,0,3,5};
  int rowc[5] = {0,1,2,0,2};
  C = value_creator(3,3,rowc,colc);
*/
  A = mread(2000000, 2000000, "col_startA.txt", "row_indA.txt");
  B = mread(2000000, 2000000, "col_startB.txt", "row_indB.txt");
  C = mread(2000000, 2000000, "col_startC.txt", "row_indC.txt");

  struct timespec ts_start, ts_end;
  clock_gettime(CLOCK_MONOTONIC, &ts_start); // get the start

  result = bmm(A, B);
  check(C,result);


  clock_gettime(CLOCK_MONOTONIC, &ts_end); // get the finishing
  long delta_sec = ts_end.tv_nsec - ts_start.tv_nsec;
  printf("Execution time: %ld ns\n", delta_sec);
/*
  printf("\nresult->row_ind = ");
  for(int i=0;i<result->col_start[C->n];i++){
    //printf("%d", i);
    printf("%d", result->row_ind[i]);
  }

  printf("\nresult->col_start = ");
  for(int i=0;i<result->n+1;i++)
    printf("%d", result->col_start[i]);
*/

  return 0;
}
