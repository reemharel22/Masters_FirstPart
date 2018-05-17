#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <omp.h>
#include "tridFunc.h"
#define N 5001
#define X 5000
double epsilon = 0.0000001;
void setMatrixForE(double (*mat)[X][X],double , double,double [N]);
void setLambda(double E[][N],double (*D)[X],double deltaX,int j);
void copyFromSolution(double*solve,double(*mat)[X][N],int j);
double integrateX(double(*U)[N][N],int i,double);
void sendToFile(double E[X][N],double T[X][N],double,double,int);
void buildDn(int n,double chi, double E[][N],double (*D)[X],double deltaX);
void buildDKershaw( double E[][N],double (*D)[X],double deltaX,int j);
void buildDMinerbo( double E[][N],double (*D)[X],double deltaX,int j);
void buildDLP( double E[][N],double (*D)[X],double deltaX,int j);
double Em[X][N],Tm[X][N],matrix[X][X];
void PredictorCorrectorFL(int times,int i,void(*)(),double,double);
double L[X],U[X],mainD[X],solve[X],D[X];
double sigmaA = 1.0;
double sigmaS = 0;
int main(void) {
  int c,k,p,h,alpha = 1,i=0,j=0;
  double lambdaE,lambdaT,tmp,x0,t0;
  double a,b,d,Src,sigmaA = 1.0;
  double sigmaS = 1.0 - sigmaA;
  double deltaX,deltaT;
  double avg = 0.0;
  FILE*fp;
  void (*funcptr) (double[][N],double(*)[X],double,int);
  //solve contains the intiail value, after solveTraigonal it will contian t+1 solution
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR Temperture ie j const, i not. you move in X axis
  alpha = 0;
  printf("Enter type of Flux limiter:\n0-Kershaw\n1-Levemore Pomraning\n2-Minerbo\n");
  scanf(" %d",&p);
  //boundry of the source
  x0 = 0.5;
  t0 = 10.0;
  c = 1.0; //speed of light
  deltaX = (double)0.01;
  deltaT = (double)0.01;
  lambdaE = (c * c * deltaT) / (deltaX*deltaX);
  lambdaT = c * deltaT * sigmaA;
  Src = 1;
  //setting up the matrices
  for (i = 0; i < N; i++) {
    //the initial function for solve aka u(x,t-1)
    Src = 1.0;
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
     for ( j = 0; j < X; j++) {
         matrix[i][j] = 0.0;
     }
  }
  //sets the matrix - Au(x,t) = u(x,t-1)
  if (p == 0) {
      funcptr = &buildDKershaw;
  } else if (p == 1) {
      funcptr = &buildDLP;
  } else if (p == 2) {
      funcptr = &buildDMinerbo;
  }

  for (i = 1; i < N; i++) {
      //double copySolution[X];
    //  for ( j = 0; j < X; j++) {
    //      copySolution[j] = solve[j]; //copySolution contains E(n),F(n)
    //  }
     // PredictorCorrectorFL(1,i,funcptr,deltaX,deltaT);
      (*funcptr)(Em,&D,deltaX,i-1);
      setMatrixForE(&matrix,deltaX,deltaT,D);
      //constructs the upper,lower diagonals.
      constructLUD(X,&L,&U,&mainD,matrix);
    //if we want to change the deltaT need TODO:
    //update lambdaE & lambdaT and call setMatrixforE & constructLUD again
    solveTriagonal(X,&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to Em.
    copyFromSolution(solve,&Em,i);
    for ( j = 0; j < X; j++) {
        //printf("%d",omp_get_num_threads());
      //this is where we calculate the next Temperture !
        Tm[j][i] = ((Tm[j][i-1] / deltaT) + sigmaA*Em[j][i]) / (sigmaA + (1.0/deltaT) );
         Src = 1;
         if (j*deltaX >= x0 || i*deltaT >= t0)
         {
             Src = 0;
         }
          solve[j] += lambdaT*Tm[j][i] + Src*deltaT; //this is for the next step
      }
    }
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
    sendToFile(Em,Tm,deltaX,deltaT,p);
    return 0;
}

double integrateX(double(*U)[N][N],int j,double deltaX){
  double sum = 0;
  int i;
  for ( i = 1; i < N; i++) {
      sum += (deltaX * ((*U)[i-1][j] + (*U)[i][j]))/2;
  }
  return sum;
}

//copies from solve to the matrix
void copyFromSolution(double*solve,double(*mat)[X][N],int j) {
  int i = 0;
  for ( i = 0; i < X; i++) {
    (*mat)[i][j] = solve[i];
  }
}

void setLambda(double E[][N],double (*D)[X],double deltaX,int j) {
    int i;
    double weff;
    double a,b;
    double R;
    for ( i = 0; i < X-1; i++) {
        if (E[i][j] != 0) {
            weff = (sigmaA*Tm[i][j] + sigmaS*E[i][j])/E[i][j];
            if (Tm[i][j] == 0.0) {
                weff = epsilon;
            }
            R = fabs((1.0 / (weff * E[i][j])) * ( ( -E[i][j] + E[i+1][j] ) / ( deltaX))) ;
        } else {
            weff = (sigmaA*Tm[i][j] + sigmaS*E[i][j])/epsilon;
            if (Tm[i][j] == 0.0) {
                weff = 1.0;
            }
            R = fabs((1.0 /(weff *(epsilon) + E[i][j])) * ( ( E[i][j] - E[i+1][j] ) / ( deltaX)) );
        }
        (*D)[i] = R;
    }
    (*D)[X-1] = (*D)[X-2];
}

void buildDKershaw(double E[][N],double (*D)[X],double deltaX,int j) {
    int i;
    double lambda1;
    double R;
    setLambda(E,D,deltaX,j);
    for (i = 0; i < X; i++) {
        double weff = (sigmaA*(Tm[i][j] + 1e-15) )/(E[i][j] + 1e-15);
        R = (*D)[i];
        lambda1 = 2.0 / ( 3.0 + sqrt(9.0 + 4.0* pow(R, 2.0) ) );
        (*D)[i] = lambda1 / weff;
    }
}

void buildDMinerbo(double E[][N],double (*D)[X],double deltaX,int j) {
    int i;
    setLambda(E,D,deltaX,j);
    for (i = 0; i < X; i++) {
        double weff = (sigmaA*(Tm[i][j] + 1e-15) )/(E[i][j] + 1e-15);
            if ((*D)[i] <= 0.9) {
                if ((*D)[i] != 0.0 ) {
                    (*D)[i] = ( 5.0 * ( sqrt(1+0.8*pow((*D)[i],2.0)) - 1.0)) / (6.0* pow((*D)[i],2.0));
                } else {
                    (*D)[i] = ( 5.0 * ( sqrt(1+0.8*pow((*D)[i] + epsilon,2.0)) - 1.0)) / (6.0* pow((*D)[i] + epsilon,2.0));
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

void buildDLP(double E[][N],double (*D)[X],double deltaX,int j) {
    int i;
    setLambda(E,D,deltaX,j);
    for (i = 0; i < X; i++) {
        double weff = (sigmaA*(Tm[i][j] + 1e-15) )/(E[i][j] + 1e-15);
        if ((*D)[i] != 0) {
            (*D)[i] = (1.0/(*D)[i]) * ( 1.0/tanh((*D)[i]) - (1.0/(*D)[i]) );
        }
        else {
            (*D)[i] = (1.0/((*D)[i] + epsilon)) * ( 1.0/(tanh((*D)[i] + epsilon) ) - (1.0/((*D)[i] + epsilon)) );
        }
        (*D)[i] = (*D)[i] / weff;
    }
}

//sets the triagonal matrix for E
void setMatrixForE(double (*mat)[X][X],double deltaX,double deltaT,double D[X]) {
  int i = 0;
  int j = 0;
  double lambda = deltaT / (deltaX*deltaX);
 /*for (i = 0; i < X; i++) {
        (*mat)[i][i] = 1 + deltaT + lambda*(D[i] + D[i-1]) ;
        if (i != 0 ) {
            (*mat)[i][i-1] = -lambda * D[i-1];
        }
        if (i != X-1) {
            (*mat)[i][i+1] = -lambda * D[i];
        }
    }*/
    //(*mat)[0][0] = 1 + deltaT + lambda*D[0];
    for (i = 0; i < X; i++) {
           (*mat)[i][i] = 1 + deltaT + lambda*(D[i+1]+2*D[i] + D[i-1])/(2.0);
           if (i != 0 ) {
               (*mat)[i][i-1] = -lambda * (D[i-1]+D[i])/2.0;
           }
           if (i != X-1) {
               (*mat)[i][i+1] = -lambda * (D[i] + D[i+1])/2.0;
           }
       }
       (*mat)[0][0] = 1 + deltaT + lambda*(D[0] + D[1])/2.0;
}

void sendToFile(double Em[X][N],double Tm[X][N],double deltaX,double deltaT,int p) {
    int i = 0,j;
    FILE*fp;
    if (p == 1)
        fp = fopen("../data/SuOlsonFluxLimitersDataLP.txt","w");
    else if (!p) {
        fp = fopen("../data/SuOlsonFluxLimitersDataKershaw.txt","w");
    }
    if (p == 2) {
        fp = fopen("../data/SuOlsonFluxLimitersDataMinerbo.txt","w");
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
      for ( i = 0; i < X; i++) {
            fprintf(fp,"%f ",Em[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
}

void PredictorCorrectorFL(int times,int i,void(*func)(),double deltaX,double deltaT) {
    int j,k;
    double copySolution[X];
    for ( j = 0; j < X; j++) {
        copySolution[j] = solve[j];//copySolution contains E(n),F(n)
    }

    (*func)(Em,&D,deltaX,i-1);
    setMatrixForE(&matrix,deltaX,deltaT,D);
    constructLUD(X,&L,&U,&mainD,matrix);
    solveTriagonal(X,&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to E.
    copyFromSolution(solve,&Em,i);
    //note, when we solve the real E(n+1) we need copySolution

copyFromSolution(solve,&Em,i);
}
