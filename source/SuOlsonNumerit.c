#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <omp.h>
#define N 5000
#define X 5000
void solveTriagonal(float(*solve)[N],float L[N],float U[N],float mainD[N]);

void constructLUD(float (*L)[N],float (*U)[N],float (*mainD)[N],float mat[N][N]);
void setMatrixForE(float (*mat)[X][X],float , float);
void copyFromSolution(float*solve,float(*mat)[X][N],int j);
void copyToSolution(float(*solve)[N],float mat[N][N],int j);
float integrateX(float(*U)[N][N],int i,float);
void sendToFile(float E[N][N],float T[N][N],float,float);
//void mallocMatrices(float**)
float Em[X][N],Tm[X][N],matrix[X][X],L[X],U[X],mainD[X],solve[X];

int main(void) {
  int c,k,p,h,i=0,j=0;
  float lambdaE,lambdaT,tmp,x0,t0,sigma,alpha;
  //float Em[N][N],Tm[N][N],matrix[N][N],L[N],U[N],mainD[N],solve[N];
  float a,b,d,Src;
  float deltaX,deltaT;
  clock_t start1,end1;
  float avg = 0.0;
  FILE*fp;
  //solve contains the intiail value, after solveTraigonal it will contian t+1 solution
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR Temperture ie j const, i not. you move in X axis
  sigma = 1.0;
  alpha = 0;
  //boundry of the source
  x0 = 0.5;
  t0 = 10.0;
  c = 1.0; //speed of light
  deltaX = (float)0.01;
  deltaT = (float)0.01;
  lambdaE = (c * c * deltaT) / (3.0*deltaX*deltaX);
  lambdaT = c * deltaT * sigma;
  Src = 1;
  //setting up the matrices
  for (i = 0; i < N; i++) {
    //the initial function for solve aka u(x,t-1)
    Src=1;
    if (i*deltaX >= x0)
    {
        Src = 0;
    }
    solve[i] = Src*deltaT;
    for ( j = 0; j < X; j++) {
        Tm[j][i] = Em[j][i] = 0.0;
    }
  }
  for ( i = 0; i < X; i++) {
     for ( j = 0; j < N; j++) {
         matrix[i][j] = 0.0;
     }
  }

  //sets the matrix - Au(x,t) = u(x,t-1)
  setMatrixForE(&matrix,lambdaE,lambdaT);
  //constructs the upper,lower diagonals.
  constructLUD(&L,&U,&mainD,matrix);
  for (i = 1; i < N; i++) {
    //if we want to change the deltaT need TODO:
    //update lambdaE & lambdaT and call setMatrixforE & constructLUD again
    solveTriagonal(&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to Em.
    copyFromSolution(solve,&Em,i);
    constructLUD(&L,&U,&mainD,matrix);
    for ( j = 0; j < X; j++) {
      //this is where we calculate the next Temperture !
        Tm[j][i] = ((Tm[j][i-1] / deltaT) + sigma*Em[j][i]) / (sigma + (1.0/deltaT) );
         Src = 1.0;
         if (j*deltaX >= x0 || i*deltaT >= t0)
         {
             Src = 0;
         }
          solve[j] += lambdaT*Tm[j][i] + Src*deltaT; //this is for the next step
    }
    }
    //printf("time took for the Tempreture calculation = %f\n",avg/N);
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < X; j++) {
                if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162|| i == 10000) {
                    if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
                        printf("%f\t",Em[j][i]);
                }
        }
        if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162|| i == 10000)
            printf("\n");
    }
    sendToFile(Em,Tm,deltaX,deltaT);
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

//copies from solve to the matrix
void copyFromSolution(float*solve,float(*mat)[X][N],int j) {
  int i = 0;
  for ( i = 0; i < X; i++) {
    (*mat)[i][j] = solve[i];
  }
}

//copies from the matrix to the solve
void copyToSolution(float(*solve)[N],float mat[N][N],int j){
  int i;
  for ( i = 0; i < N; i++) {
    (*solve)[i] = mat[i][j];
  }
}

//sets the triagonal matrix for E
void setMatrixForE(float (*mat)[X][X],float lambda,float lambda2) {
  int i = 0;
  for (i = 0; i < X; i++) {
        (*mat)[i][i] = 1 + 2*lambda + lambda2;
        if (i != 0 ) {
          (*mat)[i][i-1] = -lambda;
        }
        if (i != N-1) {
          (*mat)[i][i+1] = -lambda;
        }
  }
  (*mat)[0][0] = 1 + lambda + lambda2;
}


void sendToFile(float Em[X][N],float Tm[X][N],float deltaX,float deltaT) {
    int i = 0,j;
    FILE*fp;
    fp = fopen("../data/SuOlsonData.txt","w");
    for ( i = 0; i < X; i++) {
            fprintf(fp, "%f ",deltaX*(i) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%f ",deltaT*(i) );
    }
    fprintf(fp, "\n");
    for ( j = 0; j < N; j++) {
      for ( i = 0; i< X; i++) {
            fprintf(fp,"%f ",Em[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
}


void solveTriagonal(float(*solve)[N],float L[N],float U[N],float mainD[N]) {
    int i;
    U[0] = U[0] / mainD[0];
    (*solve)[0] = (*solve)[0] / mainD[0];

    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (i = 1; i < N; i++) {
        const float m = 1.0 / (mainD[i] - L[i] * U[i - 1]);
        U[i] = U[i] * m;
        (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
    }

    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (i = N - 2; i >=0 ; i--)
        (*solve)[i] -= U[i] * (*solve)[i + 1];
}

void constructLUD(float (*L)[N],float (*U)[N],float (*mainD)[N],float mat[N][N]){
    int i,j;
    for (i = 0; i < N-1; i++) {
      (*U)[i] = mat[i][i+1];
    }
      (*L)[0] = 0;
    for ( i = 1; i < N; i++) {
        (*L)[i] = mat[i][i-1];
    }
    #pragma omp parallel for default(shared)
    for (i = 0; i < N; i++) {
      (*mainD)[i] = mat[i][i];
    }
}
