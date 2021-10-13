#ifndef _HELPER_H_
#define _HELPER_H_

typedef struct Matrix
{
  int m;
  int n;
  int *row_ind;
  int *col_start;

} Matrix;

typedef struct Vector
{
  int n;
  int *ptr;

} Vector;

int check(Matrix* result, Matrix* C);
Matrix* mread(int m, int n, char* filename, char* filename2);
Matrix* basic_creator(int m, int n);
Matrix* value_creator(int m, int n, int* row_ind, int* col_start);
Matrix* block_creator(Matrix* m, int p, int q, int bsize);
Matrix* block2column(Matrix* m1, Matrix* b, int p);
Matrix* column2matrix(Matrix* m1, Matrix* c);
void padding(Matrix *m, char c, int p);
int binarySearch(int *vector, int start, int end, int value);
int findlocation(int *vector, int start, int end, int value);
void binaryInsertion(Vector *vector, int value);
Matrix* bma(Matrix* C1, Matrix* C2);
Matrix* bmm(Matrix* A, Matrix* B);
Matrix* bmm_omp(Matrix* A, Matrix* B, int NumMachines);
Matrix* bmm_parallel(Matrix* A, Matrix* B);
Matrix* bmmfiltered(Matrix* A, Matrix* B, Matrix* F);
Matrix* bmmfiltered_parallel(Matrix* A, Matrix* B, Matrix* F);
Matrix* bmmfiltered_omp(Matrix* A, Matrix* B, Matrix* F, int NumMachines);
Matrix* bmmblocked(Matrix* A, Matrix* B);


#endif
