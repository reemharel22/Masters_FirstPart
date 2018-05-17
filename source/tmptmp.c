#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>
#define N 1001
#define X 1001
#define NN ((X*2) + 1)
float epsilon = 0.0000001;
void printMatrix(float mat[][N]);
void printMatrix2(float mat[NN][NN]);
void setMatrixForE(float (*mat)[NN][NN],float , float,float,float [X]);
void solveTriagonal(float(*solve)[NN],float L[NN],float U[NN],float mainD[NN]);
void buildEFKershaw(float (*EF)[X], float E[X][N],float F[X+1][N],int j);
void buildEta(float (*EF)[X],float E[X][N],float F[X+1][N],int j);
void buildEFLP(float (*EF)[X], float E[X][N],float F[X+1][N],int j);
void buildEFMinerbo(float (*EF)[X], float E[X][N],float F[X+1][N],int j);
void constructLUD(float (*L)[NN],float (*U)[NN],float (*mainD)[NN],float [NN][NN]);
void copyFromSolution(float*solve,float(*mat)[X][N],float (*m)[X+1][N],int j);
void sendToFile(float E[N][N],float T[N][N],float,float,int);
void PredictorCorrectorEF(int times,int i,void(*f)(),float,float,float);
void PredictorCorrectorSolution(int times,int i, void(*f)(),float,float,float);
void CalculateT(int i,float deltaT);
void ApplyTandSource(int i,float deltaX,float deltaT,float lambdaT);
float getZ(float eta);
void BuildZ(float(*Z)[50][500]);

float E[X][N],T[X][N],matrix[NN][NN];
float L[NN],U[NN],mainD[NN],solve[NN],F[X+1][N],EF[X],Z[50][500];
float x0 = 0.5;
float t0 = 10.0;
int main() {
  int c,k,p,h,alpha = 1,i=0,j=0;
  float sigma;
  float lambdaE,lambdaT,Tp;
  float a,b,d,Src;
  float deltaX,deltaT;
  FILE*fp;
 // void (*funcptr) (float(*)[X],float[X][N],float[X][N],int);
  // now we have a trid-matrix size 2N on 2N
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR TEperture ie j const, i not. you move in X axis
  printf("Enter type of Eddington Factor:\n0-Kershaw\n1-Levemore Pomraning\n2-Minerbo\n");
  scanf(" %d",&p);
  if (!p) {
    //  funcptr = &buildEFKershaw;
  } else if (p == 1) {
    //  funcptr = &buildEFLP;
      BuildZ(&Z);
  } else if (p == 2) {
    //  funcptr = &buildEFMinerbo;
      BuildZ(&Z);
  }
  sigma = 1;
  alpha = 1.0;
  //boundry of the source
  c = 1; //speed of light
  deltaX = (float)0.01;
  deltaT = (float)0.01;
  lambdaE = deltaX;
  lambdaT = c * deltaT * sigma;
  Src = 1;
  //setting up the matrices
  for ( i = 0; i < NN; i++) {
      if (i%2 != 0) {
          Src=1;
          if (i*deltaX >= 2*x0)
          {
              Src = 0;
          }
          solve[i] = Src*deltaT;
      } else {
          solve[i] = 0;
      }
  }
  // omp parallel for default(shared) private(i,j)
  for (i = 0; i < X; i++) {
    for ( j = 0; j < N; j++) {
        T[i][j] = E[i][j] = 0.0;
    }
  }
  for (i = 0; i < X+1; i++) {
    for ( j = 0; j < N; j++) {
        F[i][j] = 0.0;
    }
  }
  // omp parallel for default(shared) private(i,j)
  for (i = 0; i < NN; i++) {
    for ( j = 0; j < NN; j++) {
        matrix[i][j] = 0.0;
    }
  }
  for ( i = 0; i < X; i++) {
      EF[i] = (float)1.0/3.0;
  }
  //sets the matrix - Au(x,t) = u(x,t-1)

  for (i = 1; i < N; i++) {
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
    PredictorCorrectorSolution(2,i, NULL,deltaX,deltaT,lambdaT);
    }
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < X; j++) {
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000) {
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
                    printf("%f\t",E[j][i]);
            }
        }
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000)
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
void solveTriagonal(float(*solve)[NN],float L[NN],float U[NN],float mainD[NN])
{
  int i;
  U[0] = U[0] / mainD[0];
  (*solve)[0] = (*solve)[0] / mainD[0];

  /* loop from 1 to X - 1 inclusive, performing the forward sweep */
  for (i = 1; i < NN; i++) {
      const float m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
      U[i] = U[i] * m;
      (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
  }

  /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */

  for (i = NN - 2; i >=0 ; i--){
      (*solve)[i] -= U[i] * (*solve)[i + 1];
  }
}

/*
  builds the L,U,mainD arrays for the Solve triagonal
*/
void constructLUD(float (*L)[NN],float (*U)[NN],float (*mainD)[NN],float mat[NN][NN]) {
  int i,j;
  // omp parallel for
  for (i = 0; i < NN-1; i++) {
    (*U)[i] = mat[i][i+1];
  }
    (*L)[0] = 0.0;
    // omp parallel for
  for ( i = 1; i < NN; i++) {
      (*L)[i] = mat[i][i-1];
     // printf("%f",(*L)[i]);
  }
  //printf("\n\n\n\n");
  // omp parallel for
  for (i = 0; i < NN; i++) {
    (*mainD)[i] = mat[i][i];
  }
}

//copies from solve to the matrix
void copyFromSolution(float*solve,float(*E)[X][N],float(*F)[X+1][N],int j) {
  int i = 0,k = 0;
  for ( i = 0; i < NN; i+=2) {
      if (k == X) {
          return;
      }
    (*E)[k][j] = solve[i+1];
    (*F)[k][j] = solve[i];
    //printf("%f\n",(*F)[k][j]);
    k++;
  }
  (*F)[k][j] = solve[NN-1];
}

/*sets the triagonal matrix for E
lambda = deltaT, lambda2 = deltaX
*/
void setMatrixForE(float (*mat)[NN][NN],float deltaX,float deltaT,float sigmaT,float EF[X]) {
  int i = 0;
  int j = 0;
  //construct the F part
  //adding the EF, need to change the F part !
 // // omp parallel for default(shared)
  for (i = 0; i < NN; i+=2) {
        (*mat)[i][i] = 1 + sigmaT;
        if (i != 0 ) {
          (*mat)[i][i-1] = (-EF[j-1]*deltaT)/(deltaX);
        }
        if (i != NN-1) {
          (*mat)[i][i+1] = (deltaT*EF[j])/(deltaX);
      }
        j++;
  }
 (*mat)[NN-2][NN-1] = (deltaT*EF[j])/(deltaX);
  //construct the E part
  for (i = 1; i < NN; i+=2) {
    (*mat)[i][i] = 1 + sigmaT;
    if (i != 1 ) {
      (*mat)[i][i-1] = -deltaT/deltaX;
    }
    if (i != NN-1) {
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

void printMatrix2(float mat[NN][NN]) {
  int i,j;
  for (i = 0; i < NN; i++) {
    for (j = 0; j < NN; j++)
      printf("%f\t", mat[i][j]);

    printf("\n");
  }
}

void sendToFile(float E[N][N],float T[N][N],float deltaX,float deltaT,int p) {
    int i = 0,j;
    FILE*fp;
    if (p == 0) {
        fp = fopen("../data/SuOlsonEddingtonFactorKershaw.txt","w");
    } else if (p == 1) {
        fp = fopen("../data/SuOlsonEddingtonFactorLP.txt","w");
    } else if (p == 2) {
        fp = fopen("../data/SuOlsonEddingtonFactorMinerbo.txt","w");
    }
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

void buildEta(float (*EF)[X], float E[X][N],float F[X+1][N],int j) {
    int i;
    for ( i = 0; i < X; i++) {
        float a = E[i][j];
        float b = F[i-1][j];
        if ( a < epsilon) {
            a = epsilon;
        }
        if ( i == 0 || i == 1) {
            b = 0;
        }
        (*EF)[i] = (float)fabsf(b/a);
    }
}

void buildEFKershaw(float (*EF)[X], float E[X][N],float F[X+1][N],int j) {
    int i;
    return;
    buildEta(EF,E,F,j);
    for (i = 0; i < X; i++) {
        (*EF)[i] = (float)((1.0 + 2.0*pow((*EF)[i],2.0))/3.0);
        if ((*EF)[i] > 1) {
            (*EF)[i] = 0.99;
        }
    }
}

void buildEFLP(float (*EF)[X], float E[X][N],float F[X+1][N],int j) {
    int i;
    buildEta(EF,E,F,j);
    if (1 == 1) {
        for (i = 0; i < X; i++) {
            float eta = (*EF)[i];
            if ((*EF)[i] < epsilon) {
                (*EF)[i] = (float)1/3;
                continue;
            } else {
                float zed = getZ((*EF)[i]);
                //printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
                if (zed == 0) {
                    zed = epsilon;
                }
                if ((*EF)[i] == 0) {
                    (*EF)[i] = epsilon;
                }
                (*EF)[i] = (-1.0/zed + 1.0/tanh(zed)) * (1.0/(tanh(zed)));
                if (fabsf((-1.0/zed + 1.0/tanh(zed))) - eta > 0.001) {
                    printf("god damnit\n");
                }
            }
        }
    } else {
        for (i = 0; i < X; i++) {
            float eta = (*EF)[i];
            (*EF)[i] = (float) (powf(eta,2.0) + 1.0/3.0 - powf(eta/1.55,2.6));
        }
    }
}

void buildEFMinerbo(float (*EF)[X], float E[X][N],float F[X+1][N],int j) {
    int i;
    if (1 == 2) {
        buildEta(EF,E,F,j);
        for (i = 0; i < X; i++) {
            if ((*EF)[i] < epsilon) {
                (*EF)[i] = (float)1/3;
                continue;
            } else {
                float zed = getZ((*EF)[i]);
                //printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
                if (zed == 0) {
                    zed = epsilon;
                }
                if ((*EF)[i] == 0) {
                    (*EF)[i] = epsilon;
                }
                (*EF)[i] =  1 - ( (2*(-1.0/zed + 1.0/tanh(zed)))/(zed));
            }
        }
    } else {
        for (i = 0; i < X; i++) {
            float b = (*EF)[i]/1.55;
            (*EF)[i] = (float)1/3 + (float)2*powf(b, 2.6);
        }
    }
}

void BuildZ(float (*Z)[50][500]) {
    int i,j;
    for ( i = 0; i < 50; i++) {
        for (j = 0; j < 500; j++) {
            if (i == 0 && j == 0) {
                (*Z)[i][j] = 0;
            } else {
                float z = j*0.002 + i;
                (*Z)[i][j] = -1.0/z + 1.0/tanh(z);
            }
        }
    }
}

float getZ(float eta) {
    int i,j;
    int k = 15;
    int l=0,r=49;

    while (1) {
        k = l + (r - l )/2;
        if ( k < 0) {
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
        if (eta > Z[k][499] && eta < Z[k+1][499]) {
            k++;
            break;
        }
        if (eta < Z[k][499] && eta > Z[k-1][499]) {
            k;
            break;
        }
        if (Z[k][499] > eta) {
            l = l;
            r = k -1;
        } else {
            r = r;
            l = k + 1;
        }
    }

    j = 250;
    r = 500;
    l = 0;
    while (1) {
        j = l + (r - l )/2;
        if ( j < 0) {
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
        if (eta > Z[k][j] && eta < Z[k][j+1]) {
            break;
        } if (eta < Z[k][j] && eta > Z[k][j-1]) {
            j;
            break;
        }
        if (Z[k][j] > eta  ) {
            l = l;
            r = j - 1;
        } else {
            r = r;
            l = j + 1;
        }
    }
    return j*0.002 + k;
}

/*
    this works with the following concept :
    we start solving with the F-L or E-F as if it is time n, which is incorrect
    we get out assumed value E* and F*, we calculate with these new values the
    E-F,F-L of time n+1, and we use them on the solution of n.
*/
void PredictorCorrectorEF(int times,int i,void(*funcptr)(),float deltaX,float deltaT,float lambdaT) {
    int j,k;
    float copySolution[NN];
    for ( j = 0; j < NN; j++) {
        copySolution[j] = solve[j];//copySolution contains E(n),F(n)
    }
    //we first do the basis, where we calculate E*,F*
    (*funcptr)(&EF,E,F,i-1);
    setMatrixForE(&matrix,deltaX,deltaT,lambdaT,EF);
    constructLUD(&L,&U,&mainD,matrix);
    solveTriagonal(&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to E.
    copyFromSolution(solve,&E,&F,i);
    //note, when we solve the real E(n+1) we need copySolution
    times--;
    while (times != 0) {
        for ( j = 0; j < NN; j++) {
            solve[j] = copySolution[j];//copySolution contains E(n),F(n)
        }
        (*funcptr)(&EF,E,F,i);
        setMatrixForE(&matrix,deltaX,deltaT,lambdaT,EF);
        constructLUD(&L,&U,&mainD,matrix);
        solveTriagonal(&solve,L,U,mainD);
        copyFromSolution(solve,&E,&F,i);
        times--;
    }
}

void PredictorCorrectorSolution(int times,int i, void(*f)(),float deltaX,float deltaT,float lambdaT) {
    int j,k;
    float En[X],Fn[X+1];
    //backing up the solutions.
    for ( j = 0; j < X; j++) {
        En[j] = E[j][i-1];
    }
    for ( j = 0; j < X+1; j++) {
        Fn[j] = F[j][i-1];
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

void CalculateT(int i,float deltaT) {
    int j,k;
    float sigma = 1;
    // omp parallel for default(shared)
    for ( j = 0; j < X; j++) {
      //this is where we calculate the next Temperture !
        T[j][i] = ((T[j][i-1] / deltaT) + sigma*E[j][i]) / (sigma + (1.0/deltaT) );
     }
}

void ApplyTandSource(int i,float deltaX,float deltaT,float lambdaT) {
    int j;
    int k = 0;
    float Src;
    for ( j = 0; j < NN; j++) {
        if ( j % 2 != 0) {
            Src = 1.0;
                if (k*deltaX > x0 || i*deltaT >= t0)
                {
                    Src = 0;
                }
                 solve[j] += lambdaT*T[k][i]+ Src*deltaT;
                 k++;
         }
    }
}
