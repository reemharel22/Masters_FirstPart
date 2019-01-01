#include "tridFunc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SuOlsonAll.h"

double epsilon1 = 0.00000000000000000000001;
double sigmaA = 1.0;
double sigmaS = 0.0;
double x1 = 0.5 ,t1 = 10;
static double Z[50][500];
static double Chi[50][500];

void solveTriagonal(int N,double(*solve)[N],double L[N],double U[N],double mainD[N]) {
    int i;
    U[0] = U[0] / mainD[0];
    (*solve)[0] = (*solve)[0] / mainD[0];

    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (i = 1; i < N; i++) {
        const double m = 1.0 / (mainD[i] - L[i] * U[i - 1]);
        U[i] = U[i] * m;
        (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
            //printf("%f\t",(L[i]));
    }

    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (i = N - 2; i >=0 ; i--) {
        (*solve)[i] -= U[i] * (*solve)[i + 1];

    }
}

void constructLUD(int N,double (*L)[N],double (*U)[N-1],double (*mainD)[N],double mat[N][N]){
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

void printMatrix(int N,double mat[][N]) {
    int i,j;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
        printf("%f\t", mat[i][j]);
      printf("\n");
    }
}

void setLambda(int X,int N,double E[][N],double Tm[][N],double (*D)[X],int j) {
    int i;
    double weff;
    double a,b;
    double R;
    double Src = 0;
    double deltaX = getdx();
    double deltaTT = getdt();
    //maybe change to Ei-1 and D[0] = 1/3
    for ( i = 0; i < X-1; i++) {
      //  weff = calculateWeff(i,j);
      //  R = fabs((1.0 / ((weff * E[i][j]) + 1e-15)) * ( ( -E[i][j] + E[i+1][j] ) / ( deltaX))) ;
        (*D)[i] = R / getOpacity(i,j);

    }
    (*D)[X-1] = (*D)[X-2];
}

void buildDKershaw(int X,double (*D)[X]) {
    /*int i;
    double lambda1;
    double R;
    double Src =0;
    double deltaX = getdx();
    double deltaTT = getdt();
    setLambda(X,N,E,Tm,D,j);
    for (i = 0; i < X; i++) {
        double weff = calculateWeff(i,j);
        R = (*D)[i];
        lambda1 = 2.0 / ( 3.0 + sqrt(9.0 + 4.0*pow(R, 2.0) ) );
        (*D)[i] = lambda1 / weff;
    }*/
}

void buildDMinerbo(int X,double (*D)[X]) {
   /* int i;
    double deltaX = getdx();
    double deltaTT = getdt();
    double Src;
    setLambda(X,N,E,Tm,D,j);
    for (i = 0; i < X; i++) {
        double weff = calculateWeff(i,j);
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
    }*/
}

void buildDLP(int X,double (*EF)[X]) {
   /* int i;
    double deltaX = getdx();
    double deltaTT = getdt();
    double Src;
    setLambda(X,N,E,Tm,D,j);
    for (i = 0; i < X; i++) {
        double weff = calculateWeff(i,j);
        if ((*D)[i] != 0) {
            (*D)[i] = (1.0/(*D)[i]) * ( 1.0/tanh((*D)[i]) - (1.0/(*D)[i]) );
        }
        else {
            (*D)[i] = (1.0/((*D)[i] + epsilon1)) * ( 1.0/(tanh((*D)[i] + epsilon1) ) - (1.0/((*D)[i] + epsilon1)) );
        }
        (*D)[i] = (*D)[i] / weff;
    }*/
}

void buildEta(int X,int N,double (*EF)[X], double E[X][N],double F[X+1][N],int j) {
    int i;
    for ( i = 0; i < X; i++) {
        double a = E[i][j];
        double b = F[i][j];
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

void buildEFKershaw(int X,double (*D)[X]) {
    /*int i;
    buildEta( X, N,EF,E,F,j);
    for (i = 0; i < X; i++) {
        double eta = (*EF)[i];
        (*EF)[i] = ((1.0 + 2.0*pow((*EF)[i],2.0))/3.0);
    }*/
}

void buildNoEF(int X,double (*D)[X])  {
    int i ;
    for ( i = 0; i < X; i++) {
       (*D[i]) = (double)1.0/(3.0 );
    }
}

void buildEFLP(int X,double (*D)[X]) {
   /* int i;
    buildEta( X, N,EF,E,F,j);
    if (1 == 2) {
        for (i = 0; i < X; i++) {
            double eta = (*EF)[i];
            if ((*EF)[i] < epsilon1) {
                (*EF)[i] = (double)1/3;
                continue;
            } else {
                double zed = getZ((*EF)[i]);
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
                double fff = 3;
            }
        }
    } else {
        for (i = 0; i < X; i++) {
            double eta = (*EF)[i];
            (*EF)[i] = (double) (pow(eta,2.0) + 1.0/3.0 - pow(eta/1.55,2.6));
        }
    }*/
}

void buildEFMinerbo(int X,double (*D)[X]) {
    int i;
  /*  buildEta( X, N,EF,E,F,j);
    if (1 == 1) {
        for (i = 0; i < X; i++) {
            double eta = (*EF)[i];
            if ((*EF)[i] < epsilon1) {
                (*EF)[i] = (double)1/3;
                continue;
            } else {
                double zed = getZ((*EF)[i]);
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
            double b = (*EF)[i]/1.55;
            (*EF)[i] = (double)1/3 + (double)2*pow(b, 2.6);
        }
    }*/
}

void buildDDiff(int X,double (*D)[X]) {
    int i;
    for ( i = 0; i < X; i++) {
        (*D)[i] = (double) 1.0/(3.0 * getOpacity(i, 0));
    }
}

void buildDDiffAsym(int X,double (*EF)[X]) {
    int i = 0;
    for (i = 0; i < X; i++) {
        //double weff = calculateWeff(i,j);
        //printf(" %lf ",weff);
        //(*D)[i] = 1.0/calculateB(weff);
    }
}

void buildDDiscDiffAsym(int X,double (*D)[X]) {
    int i = 0;
    for (i = 0; i < X; i++) {
       // double weff = calculateWeffNonAvg(i,j);
      //  (*D)[i] = 1.0/ (calculateB(weff) * getOpacity(i,j) );
      //  (*D)[i] = (*D)[i] / calculateMu2(weff);
    }
}

void BuildZ() {
    int i,j;
    for ( i = 0; i < 50; i++) {
        for (j = 0; j < 500; j++) {
            if (i == 0 && j == 0) {
                Z[i][j] = 0;
            } else {
                double z = j*0.002 + i;
                Z[i][j] = -1.0/z + 1.0/tanh(z);
            }
        }
    }
}

double getZ(double eta) {
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
