#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <mpi.h>
#include <omp.h>
#include "helper.h"
#define SERIAL_ELEMENT_LIMIT 11
#define NUM_BLOCKS 32
#define NUM_CHUNK 32
#define NUM_DIVIDE 64

Matrix* mread(int m, int n, char* filename, char* filename2)
{
  Matrix* mm;

  FILE* file;
  file = fopen(filename, "r");
  if(file == NULL)
  {
    printf("Could not open file %s.\n", filename);
    exit(-1);
  }

  int* col_start = (int *)malloc((n+1)*sizeof(int));

  for(int i=0;i<n+1;i++)
  {
    int got = fscanf(file, "%d,", &col_start[i]);
    if(got != 1)
    {
      printf("Error in reading.\n");
      exit(-2);
    }
  }
  fclose(file);

  printf("Everything is ok1\n");

  file = fopen(filename2, "r");
  if(file == NULL)
  {
    printf("Could not open file %s.\n", filename2);
    exit(-3);
  }

  int* row_ind = (int *)malloc(col_start[n]*sizeof(int));

  for(int i=0;i<col_start[n];i++)
  {
    int got = fscanf(file, "%d,", &row_ind[i]);
    if(got != 1)
    {
      printf("Error in reading.\n");
      exit(-4);
    }
  }
  fclose(file);

  mm = value_creator(m,n,row_ind,col_start);

  printf("Everything is ok2\n");
  return mm;
}

int check(Matrix* result, Matrix* C)
{
  if(result->m != C->m)
  {
    printf("Wrong result.\n");
    return -1;
  }
  else if(result->n != C->n)
  {
    printf("Wrong result.\n");
    return -1;
  }
  else
  {
    for(int i=0; i<C->n+1; i++)
    {
      if(C->col_start[i] != result->col_start[i])
      {
        printf("Wrong result.\n");
        return -1;
      }
    }

    for(int i=0; i<C->col_start[C->n]; i++)
    {
      if(C->row_ind[i] != result->row_ind[i])
      {
        printf("Wrong result.\n");
        return -1;
      }
    }

  }

  printf("Correct result.\n");
  return 0;
}

Matrix* basic_creator(int m, int n)
{

  Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
  matrix->m = m;
  matrix->n = n;
  matrix->row_ind = (int *)malloc(sizeof(int));
  matrix->col_start = (int *)malloc((n+1)*sizeof(int));

  return matrix;
}

Matrix* value_creator(int m, int n, int* row_ind, int* col_start)
{

  Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
  matrix->m = m;
  matrix->n = n;
  matrix->col_start = (int *)malloc((n+1)*sizeof(int));
  matrix->row_ind = (int *)malloc(col_start[n]*sizeof(int));

  for(int i=0;i<col_start[n];i++)
    matrix->row_ind[i] = row_ind[i];
  for(int i=0;i<n+1;i++)
    matrix->col_start[i] = col_start[i];

  return matrix;
}

Matrix* block_creator(Matrix* m, int p, int q, int bsize)
{
  Matrix* b = (Matrix*)malloc(sizeof(Matrix));
  b->m = bsize;
  b->n = bsize;
  b->col_start = (int *)malloc((bsize+1)*sizeof(int));
  b->row_ind = (int *)malloc(sizeof(int));
  b->col_start[0] = 0;

  for(int j=0; j<bsize; j++)
  {
    Vector *v = (Vector*)malloc(sizeof(Vector));
    v->n = 0;
    v->ptr = (int *)malloc(sizeof(int));
    // real j = j+q*bsize;
    numk = m->col_start[j+q*bsize+1] - m->col_start[j+q*bsize];

    for(int k=0; k<numk; k++)
    {
      int a = m->row_ind[m->col_start[j]+k];
      if(a>=p*size && a<(p+1)size)
        binaryInsertion(v, a-p*bsize);
    }

    b->col_start[j+1] = b->col_start[j] + v->n;
    b->row_ind = (int *)realloc(b->row_ind, (b->col_start[j+1]*sizeof(int)));
    for(int c=0; c<v->n; c++)
    {
      b->row_ind[C->col_start[j] + c] = v->ptr[c]-p*bsize;
    }

    free(v->ptr);
    free(v);

  }
  return b;
}

Matrix* block2column(Matrix* m1, Matrix* b, int p)
{
  Matrix* m = basic_creator(m1->m, m1->n);
  m->col_start[0] = 0;
  m->row_ind = (int *)realloc(m->row_ind, ((m->col_start[m->n+1]+b->col_start[b->n+1])*sizeof(int)));
  int bsize = b->n;

  for(int j=0; j<bsize; j++)
  {
    m->col_start[j+1] = m1->col_start[j+1] + b->col_start[j+1];
    numk1 = m1->col_start[j+1] - m1->col_start[j];
    numk = b->col_start[j+1] - b->col_start[j];

    for(int k=0; k<numk1; k++)
    {
      int a = m1->row_ind[m1->col_start[j]+k];
      m->row_ind[m->col_start[j]+k] = a;
    }

    for(int k=numk1; k<numk1+numk; k++)
    {
      int a = b->row_ind[b->col_start[j]+k-numk1] + p*bsize;
      m->row_ind[m->col_start[j]+k] = a;
    }
  }

  return m;
}

Matrix* column2matrix(Matrix* m1, Matrix* c)
{
  Matrix* m = value_creator(m1->m, m1->n, m1->row_ind, m1->col_start);
  m->n = m1->n + c->n;
  m->col_start = (int *)realloc(m->col_start, (m->n+1)*sizeof(int)));

  for(int i=0; i<c->n+1; i++){
    m->col_start[i+m1->n] = c->col_start[i] + m->col_start[m1->n];
  }

  m->row_ind = (int *)realloc(m->row_ind, (m->col_start[m->n])*sizeof(int)));
  for(int j=0; j<c->n; j++)
  {
    numk = b->col_start[j+1] - b->col_start[j];
    for(int k=0; k<numk; k++)
    {
      m->row_ind[m->col_start[j+m1->n]+k] = c->row_ind[c->col_start[j]+k];
    }
  }

  return m;
}

void padding(Matrix *m, char c, int p)
{
  if(c == 'c')
  {
    m->col_start = (int *)realloc(m->col_start, ((m->n+1 + p)*sizeof(int)));
    if(p>0)
    {
      for(int i=m->n+1; i<m->n+p+1; i++)
        m->col_start[i] = m->col_start[m->n];
    }
    m->n = m->n + p;
  }
  else if(c == 'r')
  {
    m->m = m->m + p;
  }
  else
  {
    printf("Doesn't know what kind of padding needed.\n");
    printf("Padding failed. Returning matrix unchanged..\n");
  }

}

int binarySearch(int *vector, int start, int end, int value)
{
  int numElem = end - start + 1;
  int i = start;

  if(numElem < 1)
    return 1;
  else if(numElem<SERIAL_ELEMENT_LIMIT)
  {
    for(int tmp=0; tmp<numElem; tmp++)
    {
      if(value == vector[start+tmp])
        return 0;
    }
    return 1;
  }
  else if(numElem%2 == 0)
    i = i + numElem/2;
  else
    i = i + (numElem+1)/2;


  if(value == vector[i])
    return 0;
  else if(value < vector[i])
  {
    binarySearch(vector, start, i-1, value);
  }
  else
  {
    binarySearch(vector, i+1, end, value);
  }
}

int findlocation(int *vector, int start, int end, int value)
{
  int numElem = end - start + 1;
  //printf("numElem = %d\n", numElem);
  int i = start;
  //printf("start = i = %d\n", start);
  //printf("end = %d\n", end);

  if(numElem < 1)
    return -1;
  else if(numElem<SERIAL_ELEMENT_LIMIT)
  {
    //printf("numElem < 11\n");
    while (i <= end) {
        //i = start + (end - start) / 2;
        //printf("1\n");
        if (value == vector[i]){
            //printf("2\n");
            return -1;
          }
        else if (value > vector[i]){
          i++;
          //printf("3\n");
        }
        else{return i;}
            //printf("4\n");
        }

  }
  else if(numElem%2 == 0)
    i = i + numElem/2;
  else
    i = i + (numElem+1)/2;


  if(value == vector[i])
    return -1;
  else if(value < vector[i])
  {
    findlocation(vector, start, i-1, value);
  }
  else
  {
    findlocation(vector, i+1, end, value);
  }
}

void binaryInsertion(Vector *vector, int value)
{
  if(vector->n==0)
  {
    vector->ptr[0] = value;
    vector->n++;
  }
  else
  {
    int loc = findlocation(vector->ptr, 0, vector->n-1, value);
    //printf("loc = %d\n", loc);
    if(loc != -1)
    {
      int i = vector->n-1;
      vector->n++;
      vector->ptr = (int *)realloc(vector->ptr, vector->n*sizeof(int));

      while (i >= loc)
      {
        //printf("last = %d\n", vector->ptr[i]);
        vector->ptr[i+1] = vector->ptr[i];
        i--;
      }
      vector->ptr[i+1] = value;
    }
  }


}

// boolean matrix addition
Matrix* bma(Matrix* C1, Matrix* C2)
{
  if((C1->n != C2->n) || (C1->m != C2->m))
  {
    printf("Error. Matrices cannot be added, not identical dimentions.\n");
    return NULL;
  }

  Matrix* C = basic_creator(C1->m, C1->n);
  C->col_start[0] = 0;

  for(int j=0; j<C->n; j++)
  {
    int num1 = C1->col_start[j+1] - C1->col_start[j];
    int num2 = C2->col_start[j+1] - C2->col_start[j];

    Vector *v = (Vector*)malloc(sizeof(Vector));
    v->n = num1;
    // v->ptr = C1->row_ind
    v->ptr = (int *)malloc(num1*sizeof(int));
    for(int t=0;t<num1;t++)
      v->ptr[t] = C1->row_ind[j];

    for(int i=0; i<num2; i++)
    {
      int p = C2->row_ind[C2->col_start[j]+i];
      binaryInsertion(v, p);
    }

    // this needs to be done in arithmetic order for every j
    C->col_start[j+1] = C->col_start[j] + v->n;
    C->row_ind = (int *)realloc(C->row_ind, (C->col_start[j+1]*sizeof(int)));

    for(int c=0; c<v->n; c++)
    {
      C->row_ind[C->col_start[j] + c] = v->ptr[c];
    }

    free(v->ptr);
    free(v);

  }

  return C;
}

Matrix* bmm(Matrix* A, Matrix* B)
{
  Matrix* C = basic_creator(A->m, B->n);
  C->col_start[0] = 0;

  for(int j=0; j<B->n; j++)
  {
    printf("Hi2\n");
    int numk = B->col_start[j+1] - B->col_start[j];

    Vector *v = (Vector*)malloc(sizeof(Vector));
    v->n = 0;
    v->ptr = (int *)malloc(sizeof(int));

    for(int kd=0; kd<numk; kd++)
    {
      printf("Hi3\n");
      // poia einai ta rows tou B sti stili j
      int k = B->row_ind[B->col_start[j]+kd];
      // to k einai i stili tou A
      // posa elements exei auti i stili tou A?
      printf("k[%d] = %d\n", kd, k);
      int nump = A->col_start[k+1] - A->col_start[k];
      printf("nump = %d\n", nump);
      // gia kathe k column toy pinaka A vriskw ta p stoixeia tou
      // poia einai auta ta elements? ena p gia kathe thesi pd
      for(int pd=0; pd<nump; pd++)
      {
        int p = A->row_ind[A->col_start[k]+pd];
        binaryInsertion(v, p);
        printf("nump = ERROR\n");
      }
    }
    //printf("numj[%d] = %d\n", j, v->n);
    C->col_start[j+1] = C->col_start[j] + v->n;
    C->row_ind = (int *)realloc(C->row_ind, (C->col_start[j+1]*sizeof(int)));
    //printf("%d\n", C->col_start[j+1]);
    for(int c=0; c<v->n; c++)
    {
      C->row_ind[C->col_start[j] + c] = v->ptr[c];
    }
    free(v->ptr);
    free(v);

  }
  return C;
}

Matrix* bmm_omp(Matrix* A, Matrix* B, int NumMachines)
{
  Matrix* C = basic_creator(A->m, B->n);
  C->col_start[0] = 0;

  int num_thrd = NUM_DIVIDE/NumMachines;
  if(num_thrd > 8)
    num_thrd = 8;

  int CHUNK = B->n/num_thrd;

  #pragma omp parallel for schedule(static, CHUNK) num_threads(num_thrd)
  for(int j=0; j<B->n; j++)
  {
    int numk = B->col_start[j+1] - B->col_start[j];
    //printf("numk = %d\n", numk);
    Vector *v = (Vector*)malloc(sizeof(Vector));
    v->n = 0;
    v->ptr = (int *)malloc(sizeof(int));

    for(int kd=0; kd<numk; kd++)
    {
      int k = B->row_ind[B->col_start[j]+kd];
      int nump = A->col_start[k+1] - A->col_start[k];

      for(int pd=0; pd<nump; pd++)
      {
        int p = A->row_ind[A->col_start[k]+pd];
        binaryInsertion(v, p);
        //printf("nump = ERROR\n");
      }
    }

    #pragma omp ordered
    {
      //printf("numj[%d] = %d\n", j, v->n);
      C->col_start[j+1] = C->col_start[j] + v->n;
      C->row_ind = (int *)realloc(C->row_ind, (C->col_start[j+1]*sizeof(int)));
      //printf("%d\n", C->col_start[j+1]);
      for(int c=0; c<v->n; c++)
      {
        C->row_ind[C->col_start[j] + c] = v->ptr[c];
      }
    }

    free(v->ptr);
    free(v);

  }
  return C;
}

Matrix* bmm_parallel(Matrix* A, Matrix* B)
{
  int SelfMID, NumMachines;
  MPI_Status mpistat[5];
  MPI_Request mpireq[5];

  MPI_Comm_rank( MPI_COMM_WORLD, &SelfMID );
  MPI_Comm_size( MPI_COMM_WORLD, &NumMachines );
  printf("NumTasks = %d\n", NumMachines);

  Matrix *Alocal, *partB, *partC;

  // padding B in order to have max 32 column blocks
  if( SelfMID == 0 )
  {
    int b = ceil((double)(B->n/NUM_CHUNK));
    int real_n = B->n;
    padding(B, 'c', b*NUM_CHUNK-real_n);
    // pad A as well to match dimentions
    padding(A, 'r', b*NUM_CHUNK-real_n);
  }

  // SHARING A
  if( SelfMID == 0 )
  {
    Alocal = value_creator(A->m, A->n, A->row_ind, A->col_start);
    // do with pthreads to unblock
    for(int t=1;t<NumMachines;t++)
    {
      MPI_Isend(&A->m,              1,            MPI_INT,t,0,MPI_COMM_WORLD,&mpireq[0]);
      MPI_Isend(&A->n,              1,            MPI_INT,t,1,MPI_COMM_WORLD,&mpireq[1]);
      MPI_Isend(&A->col_start[A->n],1,            MPI_INT,t,2,MPI_COMM_WORLD,&mpireq[2]);
      MPI_Isend(A->row_ind,   A->col_start[A->n], MPI_INT,t,3,MPI_COMM_WORLD,&mpireq[3]);
      MPI_Isend(A->col_start, A->n+1,             MPI_INT,t,4,MPI_COMM_WORLD,&mpireq[4]);
    }
  }
  else
  {
    int z;
    Alocal = (Matrix*)malloc(sizeof(Matrix));

    MPI_Recv(&Alocal->m, 1, MPI_INT,0,0,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&Alocal->n, 1, MPI_INT,0,1,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,         1, MPI_INT,0,2,MPI_COMM_WORLD,&mpistat[2]);

    Alocal->col_start = (int *)malloc((Alocal->n+1)*sizeof(int));
    Alocal->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(Alocal->row_ind,   z,           MPI_INT,0,3,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(Alocal->col_start, Alocal->n+1, MPI_INT,0,4,MPI_COMM_WORLD,&mpistat[4]);
  }

  // DISTRIBUTION B RECHECK
  if( SelfMID == 0 )
  {
    int pb = B->n/NumMachines;
    partB = value_creator(B->m, pb, B->row_ind, B->col_start);
    for(int t=1;t<NumMachines;t++)
    {
      int z = B->col_start[(t+1)*pb] - B->col_start[t*pb];
      //MPI_Isend(&X[t*num_points],num_points,MPI_DOUBLE,t,0,MPI_COMM_WORLD,&mpireq);
      MPI_Isend(&B->m,                             1, MPI_INT,t,5,MPI_COMM_WORLD,&mpireq[0]);
      MPI_Isend(&pb,                               1, MPI_INT,t,6,MPI_COMM_WORLD,&mpireq[1]);
      MPI_Isend(&z,                                1, MPI_INT,t,7,MPI_COMM_WORLD,&mpireq[2]);
      MPI_Isend(&B->row_ind[B->col_start[t*pb]],   z, MPI_INT,t,8,MPI_COMM_WORLD,&mpireq[3]);
      MPI_Isend(&B->col_start[t*pb],            pb+1, MPI_INT,t,9,MPI_COMM_WORLD,&mpireq[4]);
    }
  }
  else
  {
    //MPI_Recv(Xlocal,num_points,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&mpistat);
    partB = (Matrix*)malloc(sizeof(Matrix));
    int z;

    MPI_Recv(&partB->m, 1, MPI_INT,0,5,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&partB->n, 1, MPI_INT,0,6,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,        1, MPI_INT,0,7,MPI_COMM_WORLD,&mpistat[2]);

    partB->col_start = (int *)malloc((partB->n+1)*sizeof(int));
    partB->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(partB->row_ind,   z,          MPI_INT,0,8,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(partB->col_start, partB->n+1, MPI_INT,0,9,MPI_COMM_WORLD,&mpistat[4]);
  }
  // END DISTRIBUTION

  // boolean matrix multiplication
  partC = bmm_omp(Alocal, partB, NumMachines);
  //partC = bmm_mp(Alocal, partB);

  // FREEEEEEEEEEEEE partB
  free(partB->col_start);
  free(partB->row_ind);
  free(partB);

  // send back results
  if( SelfMID != 0 )
  {
    int z = partC->col_start[partC->n+1] - partC->col_start[partC->n];

    MPI_Send(&partC->m,                             1, MPI_INT,0,10,MPI_COMM_WORLD);
    MPI_Send(&partC->n,                             1, MPI_INT,0,11,MPI_COMM_WORLD);
    MPI_Send(&z,                                    1, MPI_INT,0,12,MPI_COMM_WORLD);
    MPI_Send(&B->row_ind[B->col_start[t*pb]],       z, MPI_INT,0,13,MPI_COMM_WORLD);
    MPI_Send(&B->col_start[t*pb],                pb+1, MPI_INT,0,14,MPI_COMM_WORLD);
  }
  else
  {
    int z, nz;
    // in case Irecv;
    Matrix *tmp_odd  = (Matrix*)malloc(sizeof(Matrix));
    Matrix *tmp_even = (Matrix*)malloc(sizeof(Matrix));

    // waiting to receive from t = 1
    MPI_Recv(&tmp_odd->m, 1, MPI_INT,1,10,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&tmp_odd->n, 1, MPI_INT,1,11,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,          1, MPI_INT,1,12,MPI_COMM_WORLD,&mpistat[2]);

    tmp_odd->col_start = (int *)malloc((tmp_odd->n+1)*sizeof(int));
    tmp_odd->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(tmp_odd->row_ind,   z,            MPI_INT,1,13,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(tmp_odd->col_start, tmp_odd->n+1, MPI_INT,1,14,MPI_COMM_WORLD,&mpistat[4]);

    for(int t=1;t<NumMachines-1;t++)
    {
      if(t%2 != 0)
      {
        if(t != 1)
        {
          MPI_Wait( &mpireq[3], &mpistat[3] );
          MPI_Wait( &mpireq[4], &mpistat[4] );
        }

        int nt = t+1;
        MPI_Recv(&tmp_even->m, 1, MPI_INT,nt,10,MPI_COMM_WORLD,&mpistat[0]);
        MPI_Recv(&tmp_even->n, 1, MPI_INT,nt,11,MPI_COMM_WORLD,&mpistat[1]);
        MPI_Recv(&nz,           1, MPI_INT,nt,12,MPI_COMM_WORLD,&mpistat[2]);

        tmp_even->col_start = (int *)malloc((tmp_even->n+1)*sizeof(int));
        tmp_even->row_ind = (int *)malloc(nz*sizeof(int));

        MPI_Irecv(tmp_even->row_ind,   nz,             MPI_INT,nt,13,MPI_COMM_WORLD,&mpireq[3]);
        MPI_Irecv(tmp_even->col_start, tmp_even->n+1, MPI_INT,nt,14,MPI_COMM_WORLD,&mpireq[4]);

        partC = column2matrix(partC, tmp_odd);
        free(tmp_odd->col_start);
        free(tmp_odd->row_ind);
        free(tmp_odd);
      }
      else
      {
        MPI_Wait( &mpireq[3], &mpistat[3] );
        MPI_Wait( &mpireq[4], &mpistat[4] );

        int nt = t+1;
        MPI_Recv(&tmp_odd->m, 1, MPI_INT,nt,10,MPI_COMM_WORLD,&mpistat[0]);
        MPI_Recv(&tmp_odd->n, 1, MPI_INT,nt,11,MPI_COMM_WORLD,&mpistat[1]);
        MPI_Recv(&z,        1, MPI_INT,nt,12,MPI_COMM_WORLD,&mpistat[2]);

        tmp_odd->col_start = (int *)malloc((tmp_odd->n+1)*sizeof(int));
        tmp_odd->row_ind = (int *)malloc(nz*sizeof(int));

        MPI_Irecv(tmp_odd->row_ind,   nz,          MPI_INT,nt,13,MPI_COMM_WORLD,&mpireq[3]);
        MPI_Irecv(tmp_odd->col_start, partB->n+1, MPI_INT,nt,14,MPI_COMM_WORLD,&mpireq[4]);

        partC = column2matrix(partC, tmp_even);
        free(tmp_even->col_start);
        free(tmp_even->row_ind);
        free(tmp_even);

      }

    }


    // final machine
    t = NumMachines-1;
    if(t != 1)
    {
      MPI_Wait( &mpireq[3], &mpistat[3] );
      MPI_Wait( &mpireq[4], &mpistat[4] );
    }
    partC = column2matrix(partC, tmp_odd);
    free(tmp_odd->col_start);
    free(tmp_odd->row_ind);
    free(tmp_odd);

  }

  free(Alocal->col_start);
  free(Alocal->row_ind);
  free(Alocal);

  if( SelfMID == 0 )
  {
    return partC;
  }
  else
  {
    free(partC->col_start);
    free(partC->row_ind);
    free(partC);

    return NULL;
  }

}

Matrix* bmmfiltered(Matrix* A, Matrix* B, Matrix* F)
{
  Matrix* C = basic_creator(A->m, B->n);
  C->col_start[0] = 0;

  for(int j=0; j<F->n; j++)
  {
    int numf = F->col_start[j+1] - F->col_start[j];
    int numk = B->col_start[j+1] - B->col_start[j];

    Vector *v = (Vector*)malloc(sizeof(Vector));
    v->n = 0;
    v->ptr = (int *)malloc(sizeof(int));

    if(numf!=0)
    {
      for(int kd=0; kd<numk; kd++)
      {
        int k = B->row_ind[B->col_start[j]+kd];
        //printf("OK k = %d\n", k);
        int nump  = A->col_start[k+1] - A->col_start[k];
        //int numf2 = F->col_start[k+1] - F->col_start[k];
        // gia kathe k column toy pinaka F vriskw ta p stoixeia tou
        for(int pd=0; pd<nump; pd++)
        {
          int p = A->row_ind[A->col_start[k]+pd];
          if(binarySearch(&F->row_ind[F->col_start[j]], 0, numf-1, p)==0)
          {
            printf("OK p = %d\n", p);
            binaryInsertion(v, p);
          }
        }

      }

    }
    C->col_start[j+1] = C->col_start[j] + v->n;
    C->row_ind = (int *)realloc(C->row_ind, (C->col_start[j+1])*sizeof(int));
    for(int c=0; c<v->n; c++)
    {
      C->row_ind[C->col_start[j] + c] = v->ptr[c];
    }
    // free v
    free(v->ptr);
    free(v);
  }
  return C;

}

Matrix* bmmfiltered_omp(Matrix* A, Matrix* B, Matrix* F, int NumMachines)
{
  Matrix* C = basic_creator(A->m, B->n);
  C->col_start[0] = 0;

  int num_thrd = NUM_DIVIDE/NumMachines;
  if(num_thrd > 8)
    num_thrd = 8;

  int CHUNK = B->n/num_thrd;

  #pragma omp parallel for schedule(static, CHUNK) num_threads(num_thrd)
  for(int j=0; j<F->n; j++)
  {
    int numf = F->col_start[j+1] - F->col_start[j];
    int numk = B->col_start[j+1] - B->col_start[j];

    Vector *v = (Vector*)malloc(sizeof(Vector));
    v->n = 0;
    v->ptr = (int *)malloc(sizeof(int));

    if(numf!=0)
    {
      for(int kd=0; kd<numk; kd++)
      {
        int k = B->row_ind[B->col_start[j]+kd];
        //printf("OK k = %d\n", k);
        int nump  = A->col_start[k+1] - A->col_start[k];
        //int numf2 = F->col_start[k+1] - F->col_start[k];
        // gia kathe k column toy pinaka F vriskw ta p stoixeia tou
        for(int pd=0; pd<nump; pd++)
        {
          int p = A->row_ind[A->col_start[k]+pd];
          if(binarySearch(&F->row_ind[F->col_start[j]], 0, numf-1, p)==0)
          {
            printf("OK p = %d\n", p);
            binaryInsertion(v, p);
          }
        }

      }

    }

    #pragma omp ordered
    {
      C->col_start[j+1] = C->col_start[j] + v->n;
      C->row_ind = (int *)realloc(C->row_ind, (C->col_start[j+1])*sizeof(int));
      for(int c=0; c<v->n; c++)
      {
        C->row_ind[C->col_start[j] + c] = v->ptr[c];
      }
    }

    // free v
    free(v->ptr);
    free(v);
  }
  return C;

}

Matrix* bmmfiltered_parallel(Matrix* A, Matrix* B, Matrix* F)
{
  int SelfMID, NumMachines;
  MPI_Status mpistat[5];
  MPI_Request mpireq[5];

  MPI_Comm_rank( MPI_COMM_WORLD, &SelfMID );
  MPI_Comm_size( MPI_COMM_WORLD, &NumMachines );
  printf("NumTasks = %d\n", NumMachines);

  Matrix *Alocal, *partB, *partC, *partF;

  // padding B in order to have max 32 column blocks
  if( SelfMID == 0 )
  {
    int b = ceil((double)(B->n/NUM_CHUNK));
    int real_n = B->n;
    padding(B, 'c', b*NUM_CHUNK-real_n);
    // pad A as well to match dimentions
    padding(A, 'r', b*NUM_CHUNK-real_n);
  }

  // SHARING A
  if( SelfMID == 0 )
  {
    Alocal = value_creator(A->m, A->n, A->row_ind, A->col_start);
    // do with pthreads to unblock
    for(int t=1;t<NumMachines;t++)
    {
      MPI_Isend(&A->m,              1,            MPI_INT,t,0,MPI_COMM_WORLD,&mpireq[0]);
      MPI_Isend(&A->n,              1,            MPI_INT,t,1,MPI_COMM_WORLD,&mpireq[1]);
      MPI_Isend(&A->col_start[A->n],1,            MPI_INT,t,2,MPI_COMM_WORLD,&mpireq[2]);
      MPI_Isend(A->row_ind,   A->col_start[A->n], MPI_INT,t,3,MPI_COMM_WORLD,&mpireq[3]);
      MPI_Isend(A->col_start, A->n+1,             MPI_INT,t,4,MPI_COMM_WORLD,&mpireq[4]);
    }
  }
  else
  {
    int z;
    Alocal = (Matrix*)malloc(sizeof(Matrix));

    MPI_Recv(&Alocal->m, 1, MPI_INT,0,0,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&Alocal->n, 1, MPI_INT,0,1,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,         1, MPI_INT,0,2,MPI_COMM_WORLD,&mpistat[2]);

    Alocal->col_start = (int *)malloc((Alocal->n+1)*sizeof(int));
    Alocal->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(Alocal->row_ind,   z,           MPI_INT,0,3,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(Alocal->col_start, Alocal->n+1, MPI_INT,0,4,MPI_COMM_WORLD,&mpistat[4]);
  }

  // DISTRIBUTION B RECHECK
  if( SelfMID == 0 )
  {
    int pb = B->n/NumMachines;
    partB = value_creator(B->m, pb, B->row_ind, B->col_start);
    for(int t=1;t<NumMachines;t++)
    {
      int z = B->col_start[(t+1)*pb] - B->col_start[t*pb];
      //MPI_Isend(&X[t*num_points],num_points,MPI_DOUBLE,t,0,MPI_COMM_WORLD,&mpireq);
      MPI_Isend(&B->m,                             1, MPI_INT,t,5,MPI_COMM_WORLD,&mpireq[0]);
      MPI_Isend(&pb,                               1, MPI_INT,t,6,MPI_COMM_WORLD,&mpireq[1]);
      MPI_Isend(&z,                                1, MPI_INT,t,7,MPI_COMM_WORLD,&mpireq[2]);
      MPI_Isend(&B->row_ind[B->col_start[t*pb]],   z, MPI_INT,t,8,MPI_COMM_WORLD,&mpireq[3]);
      MPI_Isend(&B->col_start[t*pb],            pb+1, MPI_INT,t,9,MPI_COMM_WORLD,&mpireq[4]);
    }
  }
  else
  {
    //MPI_Recv(Xlocal,num_points,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&mpistat);
    partB = (Matrix*)malloc(sizeof(Matrix));
    int z;

    MPI_Recv(&partB->m, 1, MPI_INT,0,5,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&partB->n, 1, MPI_INT,0,6,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,        1, MPI_INT,0,7,MPI_COMM_WORLD,&mpistat[2]);

    partB->col_start = (int *)malloc((partB->n+1)*sizeof(int));
    partB->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(partB->row_ind,   z,          MPI_INT,0,8,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(partB->col_start, partB->n+1, MPI_INT,0,9,MPI_COMM_WORLD,&mpistat[4]);
  }
  // END DISTRIBUTION

  // DISTRIBUTION F RECHECK
  if( SelfMID == 0 )
  {
    int pb = B->n/NumMachines;
    partF = value_creator(F->m, pb, F->row_ind, F->col_start);
    for(int t=1;t<NumMachines;t++)
    {
      int z = F->col_start[(t+1)*pb] - F->col_start[t*pb];
      //MPI_Isend(&X[t*num_points],num_points,MPI_DOUBLE,t,0,MPI_COMM_WORLD,&mpireq);
      MPI_Isend(&F->m,                             1, MPI_INT,t,5,MPI_COMM_WORLD,&mpireq[0]);
      MPI_Isend(&pb,                               1, MPI_INT,t,6,MPI_COMM_WORLD,&mpireq[1]);
      MPI_Isend(&z,                                1, MPI_INT,t,7,MPI_COMM_WORLD,&mpireq[2]);
      MPI_Isend(&F->row_ind[F->col_start[t*pb]],   z, MPI_INT,t,8,MPI_COMM_WORLD,&mpireq[3]);
      MPI_Isend(&F->col_start[t*pb],            pb+1, MPI_INT,t,9,MPI_COMM_WORLD,&mpireq[4]);
    }
  }
  else
  {
    //MPI_Recv(Xlocal,num_points,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&mpistat);
    partF = (Matrix*)malloc(sizeof(Matrix));
    int z;

    MPI_Recv(&partF->m, 1, MPI_INT,0,5,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&partF->n, 1, MPI_INT,0,6,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,        1, MPI_INT,0,7,MPI_COMM_WORLD,&mpistat[2]);

    partF->col_start = (int *)malloc((partF->n+1)*sizeof(int));
    partF->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(partF->row_ind,   z,          MPI_INT,0,8,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(partF->col_start, partF->n+1, MPI_INT,0,9,MPI_COMM_WORLD,&mpistat[4]);
  }
  // END DISTRIBUTION

  // boolean matrix multiplication
  partC = bmmfiltered_omp(Alocal, partB, partF, NumMachines);
  //partC = bmm_mp(Alocal, partB);
  ////////////////////

  // FREEEEEEEEEEEEE partB
  free(partB->col_start);
  free(partB->row_ind);
  free(partB);

  // FREEEEEEEEEEEEE partB
  free(partF->col_start);
  free(partF->row_ind);
  free(partF);

  // send back results
  if( SelfMID != 0 )
  {
    int z = partC->col_start[partC->n+1] - partC->col_start[partC->n];

    MPI_Send(&partC->m,                             1, MPI_INT,0,10,MPI_COMM_WORLD);
    MPI_Send(&partC->n,                             1, MPI_INT,0,11,MPI_COMM_WORLD);
    MPI_Send(&z,                                    1, MPI_INT,0,12,MPI_COMM_WORLD);
    MPI_Send(&B->row_ind[B->col_start[t*pb]],       z, MPI_INT,0,13,MPI_COMM_WORLD);
    MPI_Send(&B->col_start[t*pb],                pb+1, MPI_INT,0,14,MPI_COMM_WORLD);
  }
  else
  {
    int z, nz;
    // in case Irecv;
    Matrix *tmp_odd  = (Matrix*)malloc(sizeof(Matrix));
    Matrix *tmp_even = (Matrix*)malloc(sizeof(Matrix));

    // waiting to receive from t = 1
    MPI_Recv(&tmp_odd->m, 1, MPI_INT,1,10,MPI_COMM_WORLD,&mpistat[0]);
    MPI_Recv(&tmp_odd->n, 1, MPI_INT,1,11,MPI_COMM_WORLD,&mpistat[1]);
    MPI_Recv(&z,          1, MPI_INT,1,12,MPI_COMM_WORLD,&mpistat[2]);

    tmp_odd->col_start = (int *)malloc((tmp_odd->n+1)*sizeof(int));
    tmp_odd->row_ind = (int *)malloc(z*sizeof(int));

    MPI_Recv(tmp_odd->row_ind,   z,            MPI_INT,1,13,MPI_COMM_WORLD,&mpistat[3]);
    MPI_Recv(tmp_odd->col_start, tmp_odd->n+1, MPI_INT,1,14,MPI_COMM_WORLD,&mpistat[4]);

    for(int t=1;t<NumMachines-1;t++)
    {
      if(t%2 != 0)
      {
        if(t != 1)
        {
          MPI_Wait( &mpireq[3], &mpistat[3] );
          MPI_Wait( &mpireq[4], &mpistat[4] );
        }

        int nt = t+1;
        MPI_Recv(&tmp_even->m, 1, MPI_INT,nt,10,MPI_COMM_WORLD,&mpistat[0]);
        MPI_Recv(&tmp_even->n, 1, MPI_INT,nt,11,MPI_COMM_WORLD,&mpistat[1]);
        MPI_Recv(&nz,           1, MPI_INT,nt,12,MPI_COMM_WORLD,&mpistat[2]);

        tmp_even->col_start = (int *)malloc((tmp_even->n+1)*sizeof(int));
        tmp_even->row_ind = (int *)malloc(nz*sizeof(int));

        MPI_Irecv(tmp_even->row_ind,   nz,             MPI_INT,nt,13,MPI_COMM_WORLD,&mpireq[3]);
        MPI_Irecv(tmp_even->col_start, tmp_even->n+1, MPI_INT,nt,14,MPI_COMM_WORLD,&mpireq[4]);

        partC = column2matrix(partC, tmp_odd);
        free(tmp_odd->col_start);
        free(tmp_odd->row_ind);
        free(tmp_odd);
      }
      else
      {
        MPI_Wait( &mpireq[3], &mpistat[3] );
        MPI_Wait( &mpireq[4], &mpistat[4] );

        int nt = t+1;
        MPI_Recv(&tmp_odd->m, 1, MPI_INT,nt,10,MPI_COMM_WORLD,&mpistat[0]);
        MPI_Recv(&tmp_odd->n, 1, MPI_INT,nt,11,MPI_COMM_WORLD,&mpistat[1]);
        MPI_Recv(&z,        1, MPI_INT,nt,12,MPI_COMM_WORLD,&mpistat[2]);

        tmp_odd->col_start = (int *)malloc((tmp_odd->n+1)*sizeof(int));
        tmp_odd->row_ind = (int *)malloc(nz*sizeof(int));

        MPI_Irecv(tmp_odd->row_ind,   nz,          MPI_INT,nt,13,MPI_COMM_WORLD,&mpireq[3]);
        MPI_Irecv(tmp_odd->col_start, partB->n+1, MPI_INT,nt,14,MPI_COMM_WORLD,&mpireq[4]);

        partC = column2matrix(partC, tmp_even);
        free(tmp_even->col_start);
        free(tmp_even->row_ind);
        free(tmp_even);

      }

    }


    // final machine
    t = NumMachines-1;
    if(t != 1)
    {
      MPI_Wait( &mpireq[3], &mpistat[3] );
      MPI_Wait( &mpireq[4], &mpistat[4] );
    }
    partC = column2matrix(partC, tmp_odd);
    free(tmp_odd->col_start);
    free(tmp_odd->row_ind);
    free(tmp_odd);

  }

  free(Alocal->col_start);
  free(Alocal->row_ind);
  free(Alocal);

  if( SelfMID == 0 )
  {
    return partC;
  }
  else
  {
    free(partC->col_start);
    free(partC->row_ind);
    free(partC);

    return NULL;
  }

}

Matrix* bmmblocked(Matrix* A, Matrix* B)
{
  int max, real[3];
  real[0] = A->m;
  real[1] = A->n;
  real[2] = B->n;

  max = real[0];
  if(max<real[1])
    max = real[1];
  if(max<real[2])
    max = real[2];

  padding(A, 'r', max-A->m);
  padding(A, 'c', max-A->n);
  padding(B, 'r', max-B->m);
  padding(B, 'c', max-B->n);

  int b = ceil((double)(A->m/NUM_BLOCKS));
  max = b*NUM_BLOCKS;
  padding(A, 'r', max-A->m);
  padding(A, 'c', max-A->n);
  padding(B, 'r', max-B->m);
  padding(B, 'c', max-B->n);

  Matrix* C = basic_creator(0, 0);
  C->col_start[0] = 0;

  //
  for(int q=0; q<NUM_BLOCKS; q++)
  {
    // column q of C
    Matrix* Cc = basic_creator(0, 0);
    Cc->col_start[0] = 0;
    for(int p=0; p<NUM_BLOCKS; p++)
    {
      // find c block p q
      Matrix* Cb = basic_creator(0, 0);
      Cb->col_start[0] = 0;
      for(int s=0; s=NUM_BLOCKS; s++)
      {
        Matrix* Ab = block_creator(A, p, s, b);
        Matrix* Bb = block_creator(B, s, q, b);

        if(s==0)
        {
          Cb = bmm(Ab, Bb);
        }
        else
        {
          Cb = bma(tmp, bmm(Ab, Bb));
        }

        free(Ab->col_start);
        free(Ab->row_ind);
        free(Ab);

        free(Bb->col_start);
        free(Bb->row_ind);
        free(Bb);

      }
      Cc = block2column(Cc, Cb, p);

      free(Cb->col_start);
      free(Cb->row_ind);
      free(Cb);

    }
    C = column2matrix(C, Cc);

    free(Cc->col_start);
    free(Cc->row_ind);
    free(Cc);

  }


  // unpadding
  padding(A, 'r', real[0]-A->m);
  padding(A, 'c', real[1]-A->n);
  padding(B, 'r', real[1]-B->m);
  padding(B, 'c', real[2]-B->n);
  padding(C, 'r', real[0]-A->m);
  padding(C, 'c', real[2]-A->n);

  return C;
}
