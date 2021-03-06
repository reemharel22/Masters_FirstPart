#include "tridFunc.h"
#include <math.h>
#include <stdio.h>
#include "SuOlsonAll.h"
float epsilon1 = 0.00000000000000000000001;
float sigmaA = 1.0;
float sigmaS = 0.0;
float x1 = 0.5 ,t1 = 10;
static float Z[50][500];
static float Chi[50][500];
void solveTriagonal(int N,float(*solve)[N],float L[N],float U[N],float mainD[N]) {
    int i,j;
    U[0] = U[0] / mainD[0];
    (*solve)[0] = (*solve)[0] / mainD[0];
    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (i = 1; i < N; i++) {
        const float m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
        U[i] = (float)U[i] * m;
        (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;

  }


       //printf("\n");
    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (i = N-2; i >=0 ; i--) {
          //printf("BEFORE\tsolve: %15.15lf\tU: %15.15lf\ti is: %d\tcurrentsolve is:%15.15lf\n",(*solve)[i+1],U[i],i,(*solve)[i]);
        (*solve)[i] -= (float)(U[i] * (*solve)[i + 1]);

         //printf("AFTER\tsolve: %15.15lf\tU: %15.15lf\ti is: %d\tcurrentsolve is:%15.15lf\n",(*solve)[i+1],U[i],i,(*solve)[i]);
    }

}

void constructLUD(int N,float (*L)[N],float (*U)[N-1],float (*mainD)[N],float mat[N][N]){
    int i,j;
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

void printMatrix(int N,float mat[][N]) {
    int i,j;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
        printf("%f\t", mat[i][j]);
      printf("\n");
    }
}

void setLambda(int X,int N,float E[][N],float Tm[][N],float (*D)[X],int j) {
    int i;
    float weff;
    float a,b;
    float R;
    float Src = 0;
    float deltaX = getdx();
    float deltaTT = getdt();
    //maybe change to Ei-1 and D[0] = 1/3
    for ( i = 0; i < X-1; i++) {
        weff = calculateWeff(i,j);
        R = fabs((1.0 / ((weff * E[i][j]) + 1e-15)) * ( ( -E[i][j] + E[i+1][j] ) / ( deltaX))) ;
        (*D)[i] = R / getOpacity(i,j);

    }
    (*D)[X-1] = (*D)[X-2];
}

void buildDKershaw(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j) {
    int i;
    float lambda1;
    float R;
    float Src =0;
    float deltaX = getdx();
    float deltaTT = getdt();
    setLambda(X,N,E,Tm,D,j);
    for (i = 0; i < X; i++) {
        float weff = calculateWeff(i,j);
        R = (*D)[i];
        lambda1 = 2.0 / ( 3.0 + sqrt(9.0 + 4.0*pow(R, 2.0) ) );
        (*D)[i] = lambda1 / weff;
    }
}

void buildDMinerbo(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j) {
    int i;
    float deltaX = getdx();
    float deltaTT = getdt();
    float Src;
    setLambda(X,N,E,Tm,D,j);
    for (i = 0; i < X; i++) {
        float weff = calculateWeff(i,j);
            if ((*D)[i] <= 0.9) {
                if ((*D)[i] != 0.0 ) {
                    (*D)[i] = ( 5.0 * ( sqrt(1+0.8*pow((*D)[i],2.0)) - 1.0)) / (6.0* pow((*D)[i],2.0));
                } else {
                    (*D)[i] = ( 5.0 * ( sqrt(1+0.8*pow((*D)[i] + epsilon1,2.0)) - 1.0)) / (6.0* pow((*D)[i] + epsilon1,2.0));
                }
            }
            else if ((*D)[i] > 0.9 && (*D)[i] < 14.0) {
                (*D)[i] = -1.956*pow(10,-3) + (1.0524/(*D)[i]) -1.0412/(pow((*D)[i],1.5)) +
                0.2278/(pow((*D)[i],2)) + 0.04829/(pow((*D)[i],2.5));
            } else if ((*D)[i] >= 14) {
                (*D)[i] = (1.0/(*D)[i]) * (1.0 - (2.0/(1.0+sqrt(1.0 + 4.0*(*D)[i]))) );
            }
            (*D)[i] = (*D)[i] / weff;
    }
}

void buildDLP(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j) {
    int i;
    float deltaX = getdx();
    float deltaTT = getdt();
    float Src;
    setLambda(X,N,E,Tm,D,j);
    for (i = 0; i < X; i++) {
        float weff = calculateWeff(i,j);
        if ((*D)[i] != 0) {
            (*D)[i] = (1.0/(*D)[i]) * ( 1.0/tanh((*D)[i]) - (1.0/(*D)[i]) );
        }
        else {
            (*D)[i] = (1.0/((*D)[i] + epsilon1)) * ( 1.0/(tanh((*D)[i] + epsilon1) ) - (1.0/((*D)[i] + epsilon1)) );
        }
        (*D)[i] = (*D)[i] / weff;
    }
}

void buildEta(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],int j) {
    int i;
    for ( i = 0; i < X; i++) {
        float a = E[i][j];
        float b = F[i][j];
        if ( a == 0.0) {
            a = epsilon1;
        }
        (*EF)[i] = fabs(b)/a;
        if ((*EF)[i] > 1.0) {
            (*EF)[i] = 1.0;
        }
        if (j == 1 && (*EF)[i] != 0.0) {
        //    printf("aa %20.18f\n",b*1000);
        }
    }
}

void buildEFKershaw(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j) {
    int i;
    buildEta( X, N,EF,E,F,j);
    for (i = 0; i < X; i++) {
        float eta = (*EF)[i];
        (*EF)[i] = ((1.0 + 2.0*pow((*EF)[i],2.0))/3.0);
    }
}

void buildNoEF(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j)  {
    int i ;
    for ( i = 0; i < X; i++) {
     //   (*EF[i]) = (float)1.0/(3.0 * getOpacity(i,j));
    }
}

void buildEFLP(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j) {
    int i;
    buildEta( X, N,EF,E,F,j);
    if (1 == 2) {
        for (i = 0; i < X; i++) {
            float eta = (*EF)[i];
            if ((*EF)[i] < epsilon1) {
                (*EF)[i] = (float)1/3;
                continue;
            } else {
                float zed = getZ((*EF)[i]);
                //printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
                if (zed == 0) {
                    zed = epsilon1;
                }
                if ((*EF)[i] == 0) {
                    (*EF)[i] = epsilon1;
                }
                (*EF)[i] = (-1.0/zed + 1.0/tanh(zed)) * (1.0/(tanh(zed)));
                if (fabs((-1.0/zed + 1.0/tanh(zed))) - eta > 0.001) {
                    printf("god damnit\n");
                }
                float fff = 3;
            }
        }
    } else {
        for (i = 0; i < X; i++) {
            float eta = (*EF)[i];
            (*EF)[i] = (float) (pow(eta,2.0) + 1.0/3.0 - pow(eta/1.55,2.6));
        }
    }
}

void buildEFMinerbo(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float Tm[X][N],int j) {
    int i;
    buildEta( X, N,EF,E,F,j);
    if (1 == 1) {
        for (i = 0; i < X; i++) {
            float eta = (*EF)[i];
            if ((*EF)[i] < epsilon1) {
                (*EF)[i] = (float)1/3;
                continue;
            } else {
                float zed = getZ((*EF)[i]);
                //printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
                if (zed == 0) {
                    zed = epsilon1;
                }
                if ((*EF)[i] == 0) {
                    (*EF)[i] = epsilon1;
                }
                (*EF)[i] =  1.0 - ( (2.0*(-1.0/zed + 1.0/tanh(zed)))/(zed));
                if (fabs((-1.0/zed + 1.0/tanh(zed))) - eta > 0.001) {
                    printf("god damnit\n");
                }
            }
        }
    } else {
        for (i = 0; i < X; i++) {
            float b = (*EF)[i]/1.55;
            (*EF)[i] = (float)1/3 + (float)2*pow(b, 2.6);
        }
    }
}

void buildDDiff(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float T[X][N],int j) {
    int i;
    for ( i = 0; i < X; i++) {
        (*D)[i] = (float) 1.0/3.0;
    }
}

void buildDDiffAsym(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float T[X][N],int j) {
    int i = 0;
    for (i = 0; i < X; i++) {
        float weff = calculateWeff(i,j);
        //printf(" %lf ",weff);
        (*D)[i] = 1.0/calculateB(weff);
    }
}

void buildDDiscDiffAsym(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float T[X][N],int j) {
    int i = 0;
    for (i = 0; i < X; i++) {
        float weff = calculateWeff(i,j);
        //printf(" %lf ",weff);
        (*D)[i] = 1.0/calculateB(weff);
        (*D)[i] = (*D)[i] / calculateMu2(weff);
    }
}

void BuildZ() {
    int i,j;
    for ( i = 0; i < 50; i++) {
        for (j = 0; j < 500; j++) {
            if (i == 0 && j == 0) {
                Z[i][j] = 0;
            } else {
                float z = j*0.002 + i;
                Z[i][j] = -1.0/z + 1.0/tanh(z);
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
