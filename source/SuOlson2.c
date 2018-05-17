#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
int N = 10;
void matrixDotVector(float mat[][N],float (*vec)[N][N],int j);
void printMatrix(float mat[][N]);
void setMatrixForE(float (*mat)[N][N],float , float);
void inverseTrigMatrix(float mat[][N],float (*newmat)[N][N],float,float);
void solveTriagonal(float(*x)[N],float L[N-1],float U[N-1],float mainD[N]);
void constructLUD(float (*L)[N-1],float (*U)[N-1],float (*mainD)[N],float m[N][N]);
void copyFromSolution(float*solve,float(*mat)[N][N],int j);
void copyToSolution(float(*solve)[N],float mat[N][N],int j);
float integrateX(float(*U)[N][N],int i,float);
void setVectorsABC(float (*L)[],float (*U)[],float (*mainD)[],float,float);
//void mallocMatrices(float**)

int main() {
  int c,k,p,h,sigma = 1,alpha = 1,i=0,j=0;
  float lambdaE,lambdaT;
  float Em[N][N],Tm[N][N],matrix[N][N],a,b,d,Src;
  float L[N+1],U[N+1],mainD[N],solve[N+1],initC[N+1];
  float deltaX,deltaT;
  FILE*fp;
  //solve contains the intiail value, after solveTraigonal it will contian t+1 solution
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR Temperture ie j const, i not. you move in X axis
  sigma = 1;
  alpha = 1;
  c = 1; //speed of light
  deltaX = (float)0.1;
  deltaT = (float)0.1;
  lambdaE = (c * c * deltaT) / (3*sigma*deltaX*deltaX);
  lambdaT = c * deltaT * sigma;
  Src = c*deltaT;
  //setting up the matrices
  for (i = 0; i < N; i++) {
    //the initial function for solve aka u(x,t-1)
    initC[i] = Src;
    for ( j = 0; j < N; j++) {
        Tm[i][j] = Em[i][j] = 0.0;
        matrix[i][j] = 0.0;
  }
  initC[N] = Src;
  //sets the matrix - Au(x,t) = u(x,t-1)
//  setMatrixForE(&matrix,lambdaE,lambdaT);
}
  setVectorsABC(&L,&U,&mainD,lambdaE,deltaT);
//  printMatrix(matrix);
  //constructs the upper,lower diagonals.
//  constructLUD(&L,&U,&mainD,matrix);
  for (i = 1; i < N; i++) {
    //if we want to change the deltaT need TODO:
    //update lambdaE & lambdaT and call setMatrixforE & constructLUD again
    tridag(L,mainD,U,initC,solve,N);
    //now that we solved u(x,t+1), we will copy it to Em.
    setVectorsABC(&L,&U,&mainD,lambdaE,deltaT);
    copyFromSolution(solve,&Em,i);
  //  constructLUD(&L,&U,&mainD,matrix);
    //Em i+1 contians the solution !
    //now solving for T, not acuatl T, but V
    for ( j = 0; j < N; j++) {
      //this is where we calculate the next Temperture !
      Tm[j][i] = ((Tm[j][i-1] / deltaT) + 4*Em[j][i]) / (sigma + (1/deltaT) );
      Tm[j][i] = 0.0;
      if (i*deltaX <= 0.5 && i*deltaT <= 10) {
        //because solve triagonal first needs to have the "initial" value, we add it here
        //solve is E(x,t-1) +cdeltaT +cdeltaT*Tm^4
        //initC[j+1] += deltaT*Tm[j][i]/4 + Src;
      } else {
    //    initC[j+1] += deltaT*Tm[j][i]/4;
      }
    }
    for (j = 0; j < N; j++) {
      if (i*deltaX <= 0.5 && i*deltaT <= 10) {
        //because solve triagonal first needs to have the "initial" value, we add it here
        //solve is E(x,t-1) +cdeltaT +cdeltaT*Tm^4
        initC[j+1] = solve[j+1] + deltaT*Tm[j][i]/4 +Src;
      } else {
        initC[j+1] = solve[j+1] + deltaT*Tm[j][i]/4;
      }
    }
    }

    for ( i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        //Tm[i][j] = pow((float)Tm[i][j],(float)0.25);
      }
    }
    //    printMatrix(Tm);
    printf("\n");
    printMatrix(Em);
  d = 0;
  fp = fopen("SuOlsonData.txt","w");
    for (i = 0; i < N; i++) {
      a = integrateX(&Em,i,deltaX);
      b = integrateX(&Tm,i,deltaX);
      d = a+b;
      fprintf(fp,"%f ",d);
    }
    fprintf(fp,"\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%f ",deltaT*(i) );
    }
    fclose(fp);
    return 0;
    //copying solution to a file
    //fp = fopen("SuOlsonData.txt","w");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%f ",deltaX*(i+1) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%f ",deltaT*(i+1) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      for ( j = 0; j< N; j++) {
        fprintf(fp,"%f ",Em[j][i]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
    return 0;
}
void setVectorsABC(float (*a)[],float (*c)[],float (*b)[],float lambda,float t){
  int i,j;
  for ( i = 0; i <= N; i++) {
    (*a)[i] = 0.0;
    (*b)[i] = 0.0;
    (*c)[i] = 0.0;
  }
  for (i = 2; i < N+1; i++) {
    (*a)[i] = -lambda;
  }
  for ( i = 1; i < N; i++) {
      (*c)[i] = -lambda;
  }
  for (i = 0; i <= N; i++) {
    (*b)[i] = 1+2*lambda + t;
  }
  (*b)[1] = 1 + lambda + t;
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
void solveTriagonal(float(*solve)[N],float L[N],float U[N],float mainD[N]) {
  int i;
  U[0] = U[0] / mainD[0];
  (*solve)[0] = (*solve)[0] / mainD[0];

  /* loop from 1 to X - 1 inclusive, performing the forward sweep */
  for (i = 1; i < N; i++) {
      const float m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
      U[i] = U[i] * m;
      (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
  }

  /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
  for (i = N - 2; i >=0 ; i--)
      (*solve)[i] -= U[i] * (*solve)[i + 1];

}

/*
  builds the L,U,mainD arrays for the Solve triagonal
*/
void constructLUD(float (*L)[N],float (*U)[N-1],float (*mainD)[N],float mat[N][N])
{
  int i,j;
  (*U)[N-1] = 0;
  for (i = 0; i < N-1; i++) {
    (*U)[i] = mat[i][i+1];
  }
    (*L)[0] = 0;
  for ( i = 1; i < N; i++) {
      (*L)[i] = mat[i][i-1];
  }
  for (i = 0; i < N; i++) {
    (*mainD)[i] = mat[i][i];
  }
}

//copies from solve to the matrix
void copyFromSolution(float*solve,float(*mat)[N][N],int j) {
  int i = 0;
  for ( i = 0; i < N; i++) {
    (*mat)[i][j] = solve[i+1];
  }
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

//sets the triagonal matrix for E
void setMatrixForE(float (*mat)[N][N],float lambda,float lambda2) {
  int i = 0;
  int j = 0;
  for (i = 0; i < N; i++) {
        (*mat)[i][i] = 1 + 2*lambda + lambda2;
        if (i != 0 ) {
          (*mat)[i][i-1] = -lambda;
        }
        if (j != N-1) {
          (*mat)[i][i+1] = -lambda;
        }
  }
  (*mat)[0][0] = 1 + lambda + lambda2;
}

//inverse the trigonal matrix
void inverseTrigMatrix(float mat[][N],float (*newMat)[N][N],float lambda,float lambda2) {
  int i,j;
  float a,b;
  float theta[N+1],phi[N + 1],accum;
  theta[0] = 1;
  theta[1] = mat[0][0];
  a = 1+2*lambda + lambda2;
  b = -lambda;
  phi[N] = 1;
  phi[N-1] = a;
  for (i = 2; i <= N; i++) {
    theta[i] = a*theta[i-1] -b*b*theta[i-2];
  }

  for (i = N-2; i > 0; i--) {
      phi[i] = a*phi[i + 1]-b*b*phi[i+2];
    }
  phi[0] = mat[0][0]*phi[1] - b*b*phi[2];
  for(i=0; i < N;i++)
  {
      for(j=0; j<N;j++)
      {
        if (i == j) {
          (*newMat)[i][j] = (pow(-1, i+j)*theta[i]*phi[j+1])/theta[N];
        } else if(i < j){
          (*newMat)[i][j] = (pow(-1, i+j) * pow(b,j - i) *theta[i]
          * phi[j+1])/theta[N];
          //printf(newMat[i][j])
        } else if (i > j) {
          (*newMat)[i][j] = (pow(-1, i+j) * pow(b,i - j) *theta[j]
           * phi[i+1])/theta[N];
        }
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
