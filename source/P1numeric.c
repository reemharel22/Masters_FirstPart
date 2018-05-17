#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>
#define N 10
#define X 10
#define n (N*2 + 1)
void matrixDotVector(float mat[][N],float (*vec)[N][N],int j);
void printMatrix(float mat[][N]);
void printMatrix2(float mat[][n]);
void setMatrixForE(float (*mat)[n][n],float , float,float);
void solveTriagonal(float(*solve)[n],float L[n],float U[n],float mainD[n]);
void constructLUD(float (*L)[n],float (*U)[n],float (*mainD)[n],float [n][n]);
void copyFromSolution(float*solve,float(*mat)[X+1][N],float (*m)[X][N],int j);
void copyToSolution(float(*solve)[N],float mat[N][N],int j);
float integrateX(float(*U)[N][N],int i,float);
void sendToFile(float E[N][N],float T[N][N],float,float);
//void mallocMatrices(float**)
float E[X+1][N],T[X][N],matrix[n][n];
float L[n],U[n],mainD[n],solve[n],F[X][N];
int main() {
  int c,k,p,h,alpha = 1,i=0,j=0;
  float sigma;
  float lambdaE,lambdaT,Tp,x0,t0;
  float a,b,d,Src;
  float deltaX,deltaT;
  FILE*fp;
  //void *()
  // now we have a trid-matrix size 2N on 2N
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR TEperture ie j const, i not. you move in X axis
  sigma = 1;
  alpha = 1.0;
  //boundry of the source
  x0 = 0.5;
  t0 = 10.0;
  c = 1; //speed of light
  deltaX = (float)0.01;
  deltaT = (float)0.01;
  lambdaE = deltaX;
  lambdaT = c * deltaT * sigma;
  Src = 1;
  //setting up the matrices
  for ( i = 0; i < n; i++) {
      if (i%2 == 0) {
          Src=1;
          if (i*deltaX/2 >= x0)
          {
              Src = 0;
          }
          solve[i] = Src*deltaT;
      } else {
          solve[i] = 0;
      }
  }
  for (i = 0; i < X+1; i++) {
    for ( j = 0; j < N; j++) {
        F[i][j] = T[i][j] = E[i][j] = 0.0;
    }
  }
  for (i = 0; i < n; i++) {
    for ( j = 0; j < n; j++) {
        matrix[i][j] = 0.0;
    }
  }

  //sets the matrix - Au(x,t) = u(x,t-1)
  setMatrixForE(&matrix,deltaX,deltaT,lambdaT);
  //constructs the upper,lower diagonals.
  constructLUD(&L,&U,&mainD,matrix);
  for (i = 1; i < N; i++) {
    //if we want to change the deltaT need TODO:
    //update lambdaE & lambdaT and call seTatrixforE & constructLUD again
    solveTriagonal(&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to E.
    copyFromSolution(solve,&E,&F,i);
    constructLUD(&L,&U,&mainD,matrix);
    //#pragma omp parallel for default(shared)
    for ( j = 0; j < N; j++) {
      //this is where we calculate the next Temperture !
        T[j][i] = ((T[j][i-1] / deltaT) + sigma*E[j][i]) / (sigma + (1.0/deltaT) );
     }
     //#pragma omp parallel for default(shared) private (Src)
     for ( j = 0; j < n; j+=2) {
         Src = 1;
         if (j*deltaX/2 >= x0 || i*deltaT >= t0)
         {
             Src = 0;
         }
          solve[j] += lambdaT*T[j/2][i]+ Src*deltaT;
     }
    }

    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
        //if (i == 10/2 || i == 31/2 || i == 100/2 || i == 316/2 || i == 1000/2 || i == 3162/2  || i == 10000/2) {
        if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000){
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
                    printf("%f\t",E[j][i]);
            }
        }
            //if (i == 10/2 || i == 31/2 || i == 100/2 || i == 316/2 || i == 1000/2 || i == 3162/2  || i == 10000/2)
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000)
        printf("\n");
    }
    //printMatrix(E);
    sendToFile(E,T,deltaX,deltaT);
    return 0;
}

float integrateX(float(*U)[N][N],int j,float deltaX){
  float sum = 0;
  int i;
  for ( i = 1; i < N; i++) {
      sum += (deltaX * ((*U)[i-1][j] + (*U)[i][j]))/2;
  }
  return sum;
}

/*
  solves the equation Ax=v.
  x- initialy has the v vector, later it will have the x vector
  in our case Au(x,t+1) = u(x,t), that means x intiailly has u(x,t)
  U - the upper diagonal.
  L - the lower  diagonal
  mainD - the main Diagonal
*/
void solveTriagonal(float(*solve)[n],float L[n],float U[n],float mainD[n])
{
  int i,j;
  U[0] = U[0] / mainD[0];
  (*solve)[0] = (*solve)[0] / mainD[0];
  for (j = 0; j < n; j++) {
       // printf("%lf\t",mainD[j]);
        
    }
    // printf("\n");
  /* loop from 1 to X - 1 inclusive, performing the forward sweep */
  for (i = 1; i < n; i++) {
      const float m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
 
      U[i] = U[i] * m;
           
      (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
     // printf("%lf\t", (*solve)[i]);
      
  }

  /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
  for (i = n - 2; i >=0 ; i--){
  //    (*solve)[i] -= U[i] * (*solve)[i + 1];
       printf("BEFORE\tsolve: %15.15lf\tU: %15.15lf\ti is: %d\tcurrentsolve is:%15.15lf\n",(*solve)[i+1],U[i],i,(*solve)[i]);
        (*solve)[i] -= (U[i] * (*solve)[i + 1]);

         printf("AFTER\tsolve: %15.15lf\tU: %15.15lf\ti is: %d\tcurrentsolve is:%15.15lf\n",(*solve)[i+1],U[i],i,(*solve)[i]);
       }
    printf("\n\n");
}

/*
  builds the L,U,mainD arrays for the Solve triagonal
*/
void constructLUD(float (*L)[n],float (*U)[n],float (*mainD)[n],float mat[n][n]) {
  int i,j;
  for (i = 0; i < n-1; i++) {
    (*U)[i] = mat[i][i+1];
  }
    (*L)[0] = 0;
  for ( i = 1; i < n; i++) {
      (*L)[i] = mat[i][i-1];
  }
  for (i = 0; i < n; i++) {
    (*mainD)[i] = mat[i][i];
  }
}

//copies from solve to the matrix
void copyFromSolution(float*solve,float(*E)[X+1][N],float(*F)[X][N],int j) {
  int i = 0,k = 0;
  for ( i = 0; i < n; i+=2) {
      if (k == X) {
          break;
      }
    (*E)[k][j] = solve[i];
    (*F)[k][j] = solve[i+1];
    k++;
  }
  for (i = 0; i < n; i++) {
     // printf("%lf\t",solve[i]);
    }
   // printf("\n\n");
  (*E)[X][j] = solve[n-1];
}

//copies from the matrix to the solve
void copyToSolution(float(*solve)[N],float mat[N][N],int j){
  int i;
  for ( i = 0; i < N; i++) {
    (*solve)[i] = mat[i][j];
  }
}

/**
 *
 * CONSIDERING A TRIAGONAL MATRIX ! ! !
*/
void matrixDotVector(float mat[][N], float (*vec)[N][N],int j) {
  int i,k;
  //vec[i][j+1] refers to E(x,t+1).
  for ( i = 0; i < N; i++) {
    for ( k = 0; k < N; k++) {
  //    (*vec)[i][j+1] += mat[i][k] *vec[k][j];
    }
  }
}

/*sets the triagonal matrix for E
lambda = deltaT, lambda2 = deltaX
*/
void setMatrixForE(float (*mat)[n][n],float deltaX,float deltaT,float sigmaT) {
  int i = 0;
  int j = 0;
  for (i = 1; i < n; i+=2) {
        (*mat)[i][i] = 1+ sigmaT;
        if (i != 0 ) {
          (*mat)[i][i-1] = -deltaT/(3*deltaX);
        }
        if (i != n-1) {
          (*mat)[i][i+1] = deltaT/(3*deltaX);
      }
  }
  for (i = 0; i < n; i+=2) {
    (*mat)[i][i] = 1+ sigmaT;
    if (i != 0 ) {
      (*mat)[i][i-1] = -deltaT/deltaX;
    }
    if (i != n-1) {
      (*mat)[i][i+1] = deltaT/deltaX;
  }
  }
}


//prints the matrix
void printMatrix(float mat[][N]) {
  int i,j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%f\t", mat[i][j]);

    printf("\n");
  }
}

void printMatrix2(float mat[][n]) {
  int i,j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%f\t", mat[i][j]);

    printf("\n");
  }
}

void sendToFile(float E[N][N],float T[N][N],float deltaX,float deltaT) {
    int i = 0,j;
    FILE*fp;
    fp = fopen("../data/SuOlsonP1Data.txt","w");
    for ( i = 0; i < N; i++) {
            fprintf(fp, "%f ",deltaX*(i) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%f ",deltaT*(i) );
    }
    fprintf(fp, "\n");
    for ( j = 0; j < N; j++) {
      for ( i = 0; i < N; i++) {
            fprintf(fp,"%f ",E[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
}
