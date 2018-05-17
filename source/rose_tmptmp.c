#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>
#define N 1001
#define X 1001
#define NN ((X*2) + 1)
#include "omp.h" 
float epsilon = 0.0000001;
void printMatrix(float mat[][1001]);
void printMatrix2(float mat[1001 * 2 + 1][1001 * 2 + 1]);
void setMatrixForE(float (*mat)[1001 * 2 + 1][1001 * 2 + 1],float deltaX,float deltaT,float sigmaT,float EF[1001]);
void solveTriagonal(float (*solve)[1001 * 2 + 1],float L[1001 * 2 + 1],float U[1001 * 2 + 1],float mainD[1001 * 2 + 1]);
void buildEFKershaw(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j);
void buildEta(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j);
void buildEFLP(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j);
void buildEFMinerbo(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j);
void constructLUD(float (*L)[1001 * 2 + 1],float (*U)[1001 * 2 + 1],float (*mainD)[1001 * 2 + 1],float mat[1001 * 2 + 1][1001 * 2 + 1]);
void copyFromSolution(float *solve,float (*E)[1001][1001],float (*F)[1001 + 1][1001],int j);
void sendToFile(float E[1001][1001],float T[1001][1001],float deltaX,float deltaT,int p);
void PredictorCorrectorEF(int times,int i,void (*funcptr)(),float deltaX,float deltaT,float lambdaT);
void PredictorCorrectorSolution(int times,int i,void (*f)(),float deltaX,float deltaT,float lambdaT);
void CalculateT(int i,float deltaT);
void ApplyTandSource(int i,float deltaX,float deltaT,float lambdaT);
float getZ(float eta);
void BuildZ(float (*Z)[50][500]);
float E[1001][1001];
float T[1001][1001];
float matrix[1001 * 2 + 1][1001 * 2 + 1];
float L[1001 * 2 + 1];
float U[1001 * 2 + 1];
float mainD[1001 * 2 + 1];
float solve[1001 * 2 + 1];
float F[1001 + 1][1001];
float EF[1001];
float Z[50][500];
float x0 = 0.5;
float t0 = 10.0;

int main()
{
  int c;
  int k;
  int p;
  int h;
  int alpha = 1;
  int i = 0;
  int j = 0;
  float sigma;
  float lambdaE;
  float lambdaT;
  float Tp;
  float a;
  float b;
  float d;
  float Src;
  float deltaX;
  float deltaT;
  FILE *fp;
  void (*funcptr)(float (*)[1001], float ()[1001][1001], float ()[1001][1001], int );
// now we have a trid-matrix size 2N on 2N
//ROWS IS FOR Space ie i const, j not. you move in T axis
//COLS IS FOR TEperture ie j const, i not. you move in X axis
  printf("Enter type of Eddington Factor:\n0-Kershaw\n1-Levemore Pomraning\n2-Minerbo\n");
  scanf(" %d",&p);
  if (!p) {
    funcptr = &buildEFKershaw;
  }
   else if (p == 1) {
    funcptr = &buildEFLP;
    BuildZ(&Z);
  }
   else if (p == 2) {
    funcptr = &buildEFMinerbo;
    BuildZ(&Z);
  }
  sigma = 1;
  alpha = 1.0;
//boundry of the source
//speed of light
  c = 1;
  deltaX = ((float )0.01);
  deltaT = ((float )0.01);
  lambdaE = deltaX;
  lambdaT = c * deltaT * sigma;
  Src = 1;
//setting up the matrices
  
#pragma omp parallel for private (Src,i)
  for (i = 0; i <= 2002; i += 1) {
    if (i % 2 != 0) {
      Src = 1;
      if (i * deltaX >= 2 * x0) {
        Src = 0;
      }
      solve[i] = Src * deltaT;
    }
     else {
      solve[i] = 0;
    }
  }
// omp parallel for default(shared) private(i,j)
  
#pragma omp parallel for private (i,j)
  for (i = 0; i <= 1000; i += 1) {
    
#pragma omp parallel for private (j)
    for (j = 0; j <= 1000; j += 1) {
      T[i][j] = E[i][j] = 0.0;
    }
  }
  
#pragma omp parallel for private (i,j)
  for (i = 0; i <= 1001; i += 1) {
    
#pragma omp parallel for private (j)
    for (j = 0; j <= 1000; j += 1) {
      F[i][j] = 0.0;
    }
  }
// omp parallel for default(shared) private(i,j)
  
#pragma omp parallel for private (i,j)
  for (i = 0; i <= 2002; i += 1) {
    
#pragma omp parallel for private (j)
    for (j = 0; j <= 2002; j += 1) {
      matrix[i][j] = 0.0;
    }
  }
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 1000; i += 1) {
    EF[i] = (((float )1.0) / 3.0);
  }
//sets the matrix - Au(x,t) = u(x,t-1)
  for (i = 1; i <= 1000; i += 1) {
//if we want to change the deltaT need TODO:
//update lambdaE & lambdaT and call seTatrixforE & constructLUD again
//PredictorCorrectorEF(3,i,funcptr,deltaX,deltaT,lambdaT);
/*(*funcptr)(&EF,E,F,i-1);
    setMatrixForE(&matrix,deltaX,deltaT,lambdaT,EF);
    //constructs the upper,lower diagonals.
    constructLUD(&L,&U,&mainD,matrix);
    solveTriagonal(&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to E.
    copyFromSolution(solve,&E,&F,i);*/
//ApplyTandSource(i,deltaT,lambdaT,deltaX);
    PredictorCorrectorSolution(2,i,funcptr,deltaX,deltaT,lambdaT);
  }
  for (i = 0; i <= 1000; i += 1) {
    for (j = 0; j <= 1000; j += 1) {
      if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162 || i == 10000) {
        if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177) 
          printf("%f\t",E[j][i]);
      }
    }
    if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162 || i == 10000) 
      printf("\n");
  }
//printMatrix(E);
  sendToFile(E,T,deltaX,deltaT,p);
  return 0;
}
/*
  solves the equation Ax=v.
  x- initialy has the v vector, later it will have the x vector
  in our case Au(x,t+1) = u(x,t), that means x intiailly has u(x,t)
  U - the upper diagonal.
  L - the lower  diagonal
  mainD - the main Diagonal
*/

void solveTriagonal(float (*solve)[1001 * 2 + 1],float L[1001 * 2 + 1],float U[1001 * 2 + 1],float mainD[1001 * 2 + 1])
{
  int i;
  U[0] = U[0] / mainD[0];
  ( *solve)[0] = ( *solve)[0] / mainD[0];
/* loop from 1 to X - 1 inclusive, performing the forward sweep */
  
#pragma omp parallel for private (i)
  for (i = 1; i <= 2002; i += 1) {
    const float m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
    U[i] = U[i] * m;
    ( *solve)[i] = (( *solve)[i] - L[i] * ( *solve)[i - 1]) * m;
  }
/* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
  
#pragma omp parallel for private (i)
  for (i = 1001 * 2 + 1 - 2; i >= 0; i += -1) {
    ( *solve)[i] -= U[i] * ( *solve)[i + 1];
  }
}
/*
  builds the L,U,mainD arrays for the Solve triagonal
*/

void constructLUD(float (*L)[1001 * 2 + 1],float (*U)[1001 * 2 + 1],float (*mainD)[1001 * 2 + 1],float mat[1001 * 2 + 1][1001 * 2 + 1])
{
  int i;
  int j;
// omp parallel for
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 2001; i += 1) {
    ( *U)[i] = mat[i][i + 1];
  }
  ( *L)[0] = 0.0;
// omp parallel for
  
#pragma omp parallel for private (i)
  for (i = 1; i <= 2002; i += 1) {
    ( *L)[i] = mat[i][i - 1];
// printf("%f",(*L)[i]);
  }
//printf("\n\n\n\n");
// omp parallel for
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 2002; i += 1) {
    ( *mainD)[i] = mat[i][i];
  }
}
//copies from solve to the matrix

void copyFromSolution(float *solve,float (*E)[1001][1001],float (*F)[1001 + 1][1001],int j)
{
  int i = 0;
  int k = 0;
  for (i = 0; i <= 2002; i += 2) {
    if (k == 1001) {
      return ;
    }
    ( *E)[k][j] = solve[i + 1];
    ( *F)[k][j] = solve[i];
//printf("%f\n",(*F)[k][j]);
    k++;
  }
  ( *F)[k][j] = solve[1001 * 2 + 1 - 1];
}
/*sets the triagonal matrix for E
lambda = deltaT, lambda2 = deltaX
*/

void setMatrixForE(float (*mat)[1001 * 2 + 1][1001 * 2 + 1],float deltaX,float deltaT,float sigmaT,float EF[1001])
{
  int i = 0;
  int j = 0;
//construct the F part
//adding the EF, need to change the F part !
// // omp parallel for default(shared)
  for (i = 0; i <= 2002; i += 2) {
    ( *mat)[i][i] = 1 + sigmaT;
    if (i != 0) {
      ( *mat)[i][i - 1] = -EF[j - 1] * deltaT / deltaX;
    }
    if (i != 1001 * 2 + 1 - 1) {
      ( *mat)[i][i + 1] = deltaT * EF[j] / deltaX;
    }
    j++;
  }
  ( *mat)[1001 * 2 + 1 - 2][1001 * 2 + 1 - 1] = deltaT * EF[j] / deltaX;
//construct the E part
  
#pragma omp parallel for private (i) firstprivate (sigmaT)
  for (i = 1; i <= 2002; i += 2) {
    ( *mat)[i][i] = 1 + sigmaT;
    if (i != 1) {
      ( *mat)[i][i - 1] = -deltaT / deltaX;
    }
    if (i != 1001 * 2 + 1 - 1) {
      ( *mat)[i][i + 1] = deltaT / deltaX;
    }
  }
}
//prints the matrix

void printMatrix(float mat[][1001])
{
  int i;
  int j;
  for (i = 0; i <= 1000; i += 1) {
    for (j = 0; j <= 1000; j += 1) {
      printf("%f\t",mat[i][j]);
    }
    printf("\n");
  }
}

void printMatrix2(float mat[1001 * 2 + 1][1001 * 2 + 1])
{
  int i;
  int j;
  for (i = 0; i <= 2002; i += 1) {
    for (j = 0; j <= 2002; j += 1) {
      printf("%f\t",mat[i][j]);
    }
    printf("\n");
  }
}

void sendToFile(float E[1001][1001],float T[1001][1001],float deltaX,float deltaT,int p)
{
  int i = 0;
  int j;
  FILE *fp;
  if (p == 0) {
    fp = fopen("../data/SuOlsonEddingtonFactorKershaw.txt","w");
  }
   else if (p == 1) {
    fp = fopen("../data/SuOlsonEddingtonFactorLP.txt","w");
  }
   else if (p == 2) {
    fp = fopen("../data/SuOlsonEddingtonFactorMinerbo.txt","w");
  }
  for (i = 0; i <= 1000; i += 1) {
    fprintf(fp,"%f ",(deltaX * i));
  }
  fprintf(fp,"\n");
  for (i = 0; i <= 1000; i += 1) {
    fprintf(fp,"%f ",(deltaT * i));
  }
  fprintf(fp,"\n");
  for (j = 0; j <= 1000; j += 1) {
    for (i = 0; i <= 1000; i += 1) {
      fprintf(fp,"%f ",E[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void buildEta(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j)
{
  int i;
  for (i = 0; i <= 1000; i += 1) {
    float a = E[i][j];
    float b = F[i - 1][j];
    if (a < epsilon) {
      a = epsilon;
    }
    if (i == 0 || i == 1) {
      b = 0;
    }
    ( *EF)[i] = ((float )(fabsf(b / a)));
  }
}

void buildEFKershaw(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j)
{
  int i;
  return ;
  buildEta(EF,E,F,j);
  for (i = 0; i <= 1000; i += 1) {
    ( *EF)[i] = ((float )((1.0 + 2.0 * pow(( *EF)[i],2.0)) / 3.0));
    if (( *EF)[i] > 1) {
      ( *EF)[i] = 0.99;
    }
  }
}

void buildEFLP(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j)
{
  int i;
  buildEta(EF,E,F,j);
  if (1 == 1) {
    for (i = 0; i <= 1000; i += 1) {
      float eta = ( *EF)[i];
      if (( *EF)[i] < epsilon) {
        ( *EF)[i] = ((float )1) / 3;
        continue; 
      }
       else {
        float zed = getZ(( *EF)[i]);
//printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
        if (zed == 0) {
          zed = epsilon;
        }
        if (( *EF)[i] == 0) {
          ( *EF)[i] = epsilon;
        }
        ( *EF)[i] = ((- 1.0 / zed + 1.0 / tanh(zed)) * (1.0 / tanh(zed)));
        if ((fabsf((- 1.0 / zed + 1.0 / tanh(zed))) - eta) > 0.001) {
          printf("god damnit\n");
        }
      }
    }
  }
   else {
    for (i = 0; i <= 1000; i += 1) {
      float eta = ( *EF)[i];
      ( *EF)[i] = ((float )((powf(eta,2.0)) + 1.0 / 3.0 - (powf((eta / 1.55),2.6))));
    }
  }
}

void buildEFMinerbo(float (*EF)[1001],float E[1001][1001],float F[1001 + 1][1001],int j)
{
  int i;
  if (1 == 2) {
    buildEta(EF,E,F,j);
    for (i = 0; i <= 1000; i += 1) {
      if (( *EF)[i] < epsilon) {
        ( *EF)[i] = ((float )1) / 3;
        continue; 
      }
       else {
        float zed = getZ(( *EF)[i]);
//printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
        if (zed == 0) {
          zed = epsilon;
        }
        if (( *EF)[i] == 0) {
          ( *EF)[i] = epsilon;
        }
        ( *EF)[i] = (1 - 2 * (- 1.0 / zed + 1.0 / tanh(zed)) / zed);
      }
    }
  }
   else {
    for (i = 0; i <= 1000; i += 1) {
      float b = (( *EF)[i] / 1.55);
      ( *EF)[i] = ((float )1) / 3 + ((float )2) * powf(b,2.6);
    }
  }
}

void BuildZ(float (*Z)[50][500])
{
  int i;
  int j;
  for (i = 0; i <= 49; i += 1) {
    for (j = 0; j <= 499; j += 1) {
      if (i == 0 && j == 0) {
        ( *Z)[i][j] = 0;
      }
       else {
        float z = (j * 0.002 + i);
        ( *Z)[i][j] = (- 1.0 / z + 1.0 / tanh(z));
      }
    }
  }
}

float getZ(float eta)
{
  int i;
  int j;
  int k = 15;
  int l = 0;
  int r = 49;
  while(1){
    k = l + (r - l) / 2;
    if (k < 0) {
      k = 0;
      break; 
    }
    if (k > 49) {
      k = 29;
      break; 
    }
    if (eta == Z[k][499]) {
      break; 
    }
    if (eta > Z[k][499] && eta < Z[k + 1][499]) {
      k++;
      break; 
    }
    if (eta < Z[k][499] && eta > Z[k - 1][499]) {
      k;
      break; 
    }
    if (Z[k][499] > eta) {
      l = l;
      r = k - 1;
    }
     else {
      r = r;
      l = k + 1;
    }
  }
  j = 250;
  r = 500;
  l = 0;
  while(1){
    j = l + (r - l) / 2;
    if (j < 0) {
      j = 0;
      break; 
    }
    if (j > 499) {
      j = 499;
      break; 
    }
    if (eta == Z[k][j]) {
      break; 
    }
    if (eta > Z[k][j] && eta < Z[k][j + 1]) {
      break; 
    }
    if (eta < Z[k][j] && eta > Z[k][j - 1]) {
      j;
      break; 
    }
    if (Z[k][j] > eta) {
      l = l;
      r = j - 1;
    }
     else {
      r = r;
      l = j + 1;
    }
  }
  return (j * 0.002 + k);
}
/*
    this works with the following concept :
    we start solving with the F-L or E-F as if it is time n, which is incorrect
    we get out assumed value E* and F*, we calculate with these new values the
    E-F,F-L of time n+1, and we use them on the solution of n.
*/

void PredictorCorrectorEF(int times,int i,void (*funcptr)(),float deltaX,float deltaT,float lambdaT)
{
  int j;
  int k;
  float copySolution[2002 + 1];
  
#pragma omp parallel for private (j)
  for (j = 0; j <= 2002; j += 1) {
//copySolution contains E(n),F(n)
    copySolution[j] = solve[j];
  }
//we first do the basis, where we calculate E*,F*
  ( *funcptr)(&EF,E,F,i - 1);
  setMatrixForE(&matrix,deltaX,deltaT,lambdaT,EF);
  constructLUD(&L,&U,&mainD,matrix);
  solveTriagonal(&solve,L,U,mainD);
//now that we solved u(x,t+1), we will copy it to E.
  copyFromSolution(solve,&E,&F,i);
//note, when we solve the real E(n+1) we need copySolution
  times--;
  while(times != 0){
    
#pragma omp parallel for private (j)
    for (j = 0; j <= 2002; j += 1) {
//copySolution contains E(n),F(n)
      solve[j] = copySolution[j];
    }
    ( *funcptr)(&EF,E,F,i);
    setMatrixForE(&matrix,deltaX,deltaT,lambdaT,EF);
    constructLUD(&L,&U,&mainD,matrix);
    solveTriagonal(&solve,L,U,mainD);
    copyFromSolution(solve,&E,&F,i);
    times--;
  }
}

void PredictorCorrectorSolution(int times,int i,void (*f)(),float deltaX,float deltaT,float lambdaT)
{
  int j;
  int k;
  float En[1001];
  float Fn[1001 + 1];
//backing up the solutions.
  
#pragma omp parallel for private (j)
  for (j = 0; j <= 1000; j += 1) {
    En[j] = E[j][i - 1];
  }
  
#pragma omp parallel for private (j)
  for (j = 0; j <= 1001; j += 1) {
    Fn[j] = F[j][i - 1];
  }
  PredictorCorrectorEF(1,i,f,deltaX,deltaT,lambdaT);
  CalculateT(i,deltaT);
  ApplyTandSource(i,deltaX,deltaT,lambdaT);
//after this step at E contains En+1 and T contains n+1 etc.
//and not to forget solve contains En+1 and Fn+1
//we apply to solve the addition of Tn+1 and the source
//now we will want to solve the system again with Tn+1 instead of Tn
//we copy the old E,F aka En and Fn to solve.
/*k = 0;
    int p = 0;
    for (j = 0; j < NN; j++) {
        if (j % 2 != 0) {
            solve[j] = En[k];
            k++;
        } else {
            solve[j] = Fn[p];
            p++;
        }
    }
    //now solve holds E[n] and F[n]. we apply our source and T[n+1].
    ApplyTandSource(i,deltaX,deltaT,lambdaT);
    //and now we solve again
    PredictorCorrectorEF(3,i,f,deltaX,deltaT,lambdaT);
    //now we have the correct E and F. or the more correct. now we will re-calculate
    //T and apply the source again.
    CalculateT(i,deltaT);
    ApplyTandSource(i,deltaX,deltaT,lambdaT);*/
}

void CalculateT(int i,float deltaT)
{
  int j;
  int k;
  float sigma = 1;
// omp parallel for default(shared)
  
#pragma omp parallel for private (j) firstprivate (i,deltaT,sigma)
  for (j = 0; j <= 1000; j += 1) {
//this is where we calculate the next Temperture !
    T[j][i] = ((T[j][i - 1] / deltaT + sigma * E[j][i]) / (sigma + 1.0 / deltaT));
  }
}

void ApplyTandSource(int i,float deltaX,float deltaT,float lambdaT)
{
  int j;
  int k = 0;
  float Src;
  
#pragma omp parallel for private (k,Src,j)
  for (j = 0; j <= 2002; j += 1) {
    if (j % 2 != 0) {
      Src = 1.0;
      if (k * deltaX > x0 || i * deltaT >= t0) {
        Src = 0;
      }
      solve[j] += lambdaT * T[k][i] + Src * deltaT;
      k++;
    }
  }
}
