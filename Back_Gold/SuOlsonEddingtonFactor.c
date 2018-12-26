#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>
#define N 3000
#define X 3000
#define NN (((X*2) + 1))
double epsilon = 0.00000000000000000000001;

void setMatrixForE(double (*mat)[NN][NN],double , double,double [X]);
void solveTriagonal(double(*solve)[NN],double L[NN],double U[NN],double mainD[NN]);
void buildEFKershaw(double (*EF)[X], double E[X][N],double F[X+1][N],int j);
void buildEta(double (*EF)[X],double E[X][N],double F[X+1][N],int j);
void buildEFLP(double (*EF)[X], double E[X][N],double F[X+1][N],int j);
void buildEFMinerbo(double (*EF)[X], double E[X][N],double F[X+1][N],int j);
void buildNoEF(double (*EF)[X], double E[X][N],double F[X+1][N],int j);
void buildABLambdaT(double (*EF)[X], double E[X][N],double F[X+1][N],int j);
void constructLUD(double (*L)[NN],double (*U)[NN],double (*mainD)[NN],double [NN][NN]);
void copyFromSolution(double*solve,double(*mat)[X][N],double (*m)[X+1][N],int j);
void sendToFile(double E[N][N],double T[N][N],double,double,int);
void PredictorCorrectorEF(int times,int i,void(*f)(),double,double);
void PredictorCorrectorSolution(int times,int i, void(*f)(),double,double);
void CalculateT(int i,double deltaT);
void ApplyTandSource(int i,double deltaX,double deltaT);
double getZ(double eta);
double getOpacity(int space,int time1);
void BuildZ();
int checkConverged(int i);

double E[X][N],T[X][N],matrix[NN][NN];
double L[NN],U[NN],mainD[NN],solve[NN],F[X+1][N],EF[X],Z[50][500];
double x0 = 0.5;
double t0 = 10.0;
double A[X];
int timestep = 0;
double B[X];
//double A = 3,B = 3;
double lambdaT;
double deltaX = 0.01;
double deltaT = 0.01;
int main() {
  int k,p,h,i=0,j=0;
  double sigma,alpha,c;
  double lambdaE,Tp;
  double a,b,d,Src;
  FILE*fp;
  void (*funcptr) (double(*)[X],double[X][N],double[X][N],int);
  // now we have a trid-matrix size 2N on 2N
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR TEperture ie j const, i not. you move in X axis
  printf("Enter type of Eddington Factor:\n0-Kershaw\n1-Levemore Pomraning\n2-Minerbo\n3-P1\n4-P1AB\n");
  scanf(" %d",&p);
  if (!p) {
      funcptr = &buildEFKershaw;
  } else if (p == 1) {
      funcptr = &buildEFLP;
      BuildZ(&Z);
  } else if (p == 2) {
      funcptr = &buildEFMinerbo;
      BuildZ(&Z);
  } else if (p == 3) {
    funcptr = &buildNoEF;
  } else if (p == 4) {
    funcptr = &buildABLambdaT;
  }
  sigma = 1.0;
  alpha = 1.0;
  //boundry of the source
  c = 1.0; //speed of light
  lambdaE = deltaX;
  lambdaT = (c * deltaT * sigma);
  Src = 0;
  //setting up the matrices
  for ( i = 0; i < NN; i++) {
      solve[i] = 0;
  }
  for ( i = 0; i < X; i++) {
      Src = 0;
      if (i*deltaX >= x0) {
          Src = 0;
      }
      solve[2*i + 1] = Src*deltaT;
  }
  solve[0] = 0.5;
  //#pragma omp parallel for default(shared)
  for (i = 0; i < X; i++) {
    for ( j = 0; j < N; j++) {
        T[i][j] = E[i][j] = pow(10,-5);
    }
  }
  for (i = 0; i < X+1; i++) {
    for ( j = 0; j < N; j++) {
        F[i][j] = 0.0;
    }
  }
  //#pragma omp parallel for default(shared) private(i,j)
  for (i = 0; i < NN; i++) {
    for ( j = 0; j < NN; j++) {
        matrix[i][j] = 0.0;
    }
  }

  for ( i = 0; i < X; i++) {
      EF[i] = (double)1.0/3.0;
  }
  for ( i = 0; i < X; i++) {
    A[i] = B[i] = 3.0;
  }
  for (i = 1; i < N; i++) {
//    if (p == 0 || p == 1 || p == 2) {
 
      PredictorCorrectorSolution(2,i, funcptr,deltaX,deltaT);
         timestep = i;
//    } else {
//        PredictorCorrectorSolution(1,i, funcptr,deltaX,deltaT);
//    }

    }
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < X; j++) {
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000) {
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
                   // printf("%f\t",E[j][i]);
                    printf("%f\t",pow(T[j][i],0.25));
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
void solveTriagonal(double(*solve)[NN],double L[NN],double U[NN],double mainD[NN]){
  int i;
  U[0] = U[0] / mainD[0];
  (*solve)[0] = (*solve)[0] / mainD[0];
    
  /* loop from 1 to X - 1 inclusive, performing the forward sweep */
  for (i = 1; i < NN; i++) {
      const double m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
      U[i] = U[i] * m;
    //     printf("%lf\t",(*solve)[i - 1]);
      (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
    
     //  printf("%lf\t", (*solve)[i]);
    }
   // printf("\n");

  /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
  for (i = NN - 2; i >=0 ; i--){
      (*solve)[i] -= U[i] * (*solve)[i + 1];
  }

}

/*
  builds the L,U,mainD arrays for the Solve triagonal
*/
void constructLUD(double (*L)[NN],double (*U)[NN],double (*mainD)[NN],double mat[NN][NN]) {
  int i,j;
 // #pragma omp parallel for
  for (i = 0; i < NN-1; i++) {
    (*U)[i] = mat[i][i+1];
  }
    (*L)[0] = 0.0;
  //    #pragma omp parallel for
  for ( i = 1; i < NN; i++) {
      (*L)[i] = mat[i][i-1];
     // printf("%f",(*L)[i]);
  }
  //printf("\n\n\n\n");
 //  #pragma omp parallel for
  for (i = 0; i < NN; i++) {
    (*mainD)[i] = mat[i][i];
  }
}

//copies from solve to the matrix
void copyFromSolution(double*solve,double(*E)[X][N],double(*F)[X+1][N],int j) {
  int i = 0,k = 0;
  for ( i = 0; i < NN; i+=2) {
      if (k == X) {
          return;
      }
    (*E)[k][j] = solve[i+1];
    (*F)[k][j] = solve[i];
    k++;
  }
  (*F)[X][j] = solve[NN-1];
  (*F)[0][j] = 0.5;
}

/*sets the triagonal matrix for E
lambda = deltaT, lambda2 = deltaX
*/
void setMatrixForE(double (*mat)[NN][NN],double deltaX,double deltaT,double EF[X]) {
  int i = 0;
  int j = 0;
  //construct the F part
  //adding the EF, need to change the F part !
  for (i = 0; i < NN; i+=2) {
        (*mat)[i][i] = 1 + ((deltaT * B[j])/A[j]);
        if ( j >= X) {
            (*mat)[i][i] = 1 + getOpacity(j,timestep)*((deltaT * B[X-1]) / A[X-1]);
        }
        if (i != 0 ) {
         (*mat)[i][i-1] = (-EF[j-1]*deltaT)/(deltaX);
        }
        if (i != NN-1 && i != 0) {
          (*mat)[i][i+1] = (deltaT*EF[j])/(deltaX);
      }
        j++;
  }
  //construct the E part
  j = 0;
  for (i = 1; i < NN; i+=2) {
    (*mat)[i][i] = 1 + getOpacity(j,timestep)*deltaT;
    if (i != 0 ) {
      (*mat)[i][i-1] = -deltaT/deltaX;
    }
      (*mat)[i][i+1] = deltaT/deltaX;
      j++;
  }
  (*mat)[0][0] = 1.0;
  (*mat)[0][1] = 0.5;
}

void sendToFile(double E[N][N],double T[N][N],double deltaX,double deltaT,int p) {
    int i = 0,j;
    FILE*fp;
    if (p == 0) {
        fp = fopen("../data/SuOlsonEddingtonFactorKershaw.txt","w");
    } else if (p == 1) {
        fp = fopen("../data/SuOlsonEddingtonFactorLP.txt","w");
    } else if (p == 2) {
        fp = fopen("../data/SuOlsonEddingtonFactorMinerbo.txt","w");
    } else if (p == 3 ) {
      fp = fopen("../data/Temp/SuOlsonP1Data.txt","w");
    } else if (p == 4) {
      fp = fopen("../data/SuOlsonP1AB.txt","w");
  } else if (p == 5) {

  }else if ( p == 6 ) {}
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
            fprintf(fp,"%f ",pow(T[i][j],0.25));
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
}

void buildEta(double (*EF)[X], double E[X][N],double F[X+1][N],int j) {
    int i;
    for ( i = 0; i < X; i++) {
        double a = E[i][j];
        double b = F[i][j];
        if ( a == 0.0) {
            a = epsilon;
        }
        (*EF)[i] = fabs(b)/a;
        if ((*EF)[i] > 1.0) {
            (*EF)[i] = 1.0;
        }
        if (j == 1 && (*EF)[i] != 0.0) {
        //    printf("aa %20.18f\n",b*1000);
        }
            if ((a > E[i-1][j] + 1e-5) && j == 316 && (i>1 && i<500))
                printf("%d\t%20.18f\t%20.18f\n",i,a,E[i-1][j]);
    }
}

void buildEFKershaw(double (*EF)[X], double E[X][N],double F[X+1][N],int j) {
    int i;
    buildEta(EF,E,F,j);
    for (i = 0; i < X; i++) {
        double eta = (*EF)[i];
        (*EF)[i] = ((1.0 + 2.0*pow((*EF)[i],2.0))/3.0);
    }
}

void buildNoEF(double (*EF)[X], double E[X][N],double F[X+1][N],int j)  {
  return;
}

void buildEFLP(double (*EF)[X], double E[X][N],double F[X+1][N],int j) {
    int i;
    buildEta(EF,E,F,j);
    if (1 == 2) {
        for (i = 0; i < X; i++) {
            double eta = (*EF)[i];
            if ((*EF)[i] < epsilon) {
                (*EF)[i] = (double)1/3;
                continue;
            } else {
                double zed = getZ((*EF)[i]);
                //printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
                if (zed == 0) {
                    zed = epsilon;
                }
                if ((*EF)[i] == 0) {
                    (*EF)[i] = epsilon;
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
    }
}

void buildEFMinerbo(double (*EF)[X], double E[X][N],double F[X+1][N],int j) {
    int i;
    if (1 == 1) {
        buildEta(EF,E,F,j);
        for (i = 0; i < X; i++) {
            double eta = (*EF)[i];
            if ((*EF)[i] < epsilon) {
                (*EF)[i] = (double)1/3;
                continue;
            } else {
                double zed = getZ((*EF)[i]);
                //printf("zed is %f\t eta is %f\n",zed,(*EF)[i]);
                if (zed == 0) {
                    zed = epsilon;
                }
                if ((*EF)[i] == 0) {
                    (*EF)[i] = epsilon;
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
    }
}

void buildABLambdaT(double (*EF)[X], double E[X][N],double F[X+1][N],int j) {
  //we update first lambdaT, then EF
  //EF = 1/A.
  //lambdaT = B/A ! ! !
  int i,k;
  double Src;
  for (i = 0; i < X; i++) {
      double weff;
      Src = 1.0;
      if (j*deltaT >= t0 || i*deltaX >= x0) {
        Src = 0.0;
      }
      double tt,ee;
      tt = T[i][j];
      ee = E[i][j];
      if (ee < epsilon) {
        ee = epsilon;
    }
        weff = (tt + Src ) / ee;
      if (0.55 <= weff && weff <= 0.65) {
        A[i] = 0.96835 - 0.437*weff;
      } else {
          double a = ( 0.247 * (0.433 + 0.421*weff -2.681*weff*weff
               - 1.82*pow(weff,3.0) + 4.9*pow(weff,4.0) -1.06*pow(weff,5.0)
                + 2.56*pow(weff,6.0) ) );
          double b =  pow(0.33 + 0.159*weff - 0.567 * pow(weff,2.0) - pow(weff,3.0)  ,2.0);
        A[i] = a/b ;
      }

      //we setup B now
      if ( 0.59 <= weff && weff <=0.61) {
        B[i] = 1.0 / (0.80054 - 0.523*weff);
      } else {
          double c = (1.0 + weff) / 0.40528473;
          double a = (0.1326495 + weff*(0.03424169 + weff*(0.1774006 - weff)));
          double b = (0.3267567 + weff*(0.1587312 - weff*(0.5665676 + weff)));
          B[i] = a*c/b;
      }
      if (A[i] < epsilon) {
          A[i] = epsilon;
      }
      (*EF[i]) = 1.0/A[i];
      if (j == 100 && i >30 && i < 200) {
          printf("%f\t%d\t%f\t%f\n",weff,i, (*EF[i]),B[i]/A[i]);
      }
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

/*
    this works with the following concept :
    we start solving with the F-L or E-F as if it is time n, which is incorrect
    we get out assumed value E* and F*, we calculate with these new values the
    E-F,F-L of time n+1, and we use them on the solution of n.
*/
void PredictorCorrectorEF(int times,int i,void(*funcptr)(),double deltaX,double deltaT) {
    int j,k=0,p=0;
    double copySolution[NN];
    for (j = 0; j < NN; j++) {
        if (j % 2 != 0) {
            copySolution[j] = E[k][i-1];
            k++;
        } else {
            copySolution[j] = F[p][i-1];
            p++;
        }
    }
    //we first do the basis, where we calculate E*,F*
    (*funcptr)(&EF,E,F,i-1);
    setMatrixForE(&matrix,deltaX,deltaT,EF);
    constructLUD(&L,&U,&mainD,matrix);
    solveTriagonal(&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to E.
    copyFromSolution(solve,&E,&F,i);
    //note, when we solve the real E(n+1) we need copySolution
    CalculateT(i,deltaT);//we calculate Tn+1
    ApplyTandSource(i,deltaX,deltaT);//we apply to solve Tn+1 and the src for the next step
    times--;
    times = 0;
    while (times) {
        times++;
        for ( j = 0; j < NN; j++) {
            solve[j] = copySolution[j];//copySolution contains E(n),F(n)
        }
        ApplyTandSource(i,deltaX,deltaT);
        (*funcptr)(&EF,E,F,i);
        setMatrixForE(&matrix,deltaX,deltaT,EF);//build a new more correct matrix
        constructLUD(&L,&U,&mainD,matrix);
        solveTriagonal(&solve,L,U,mainD);//solve En+1 and Fn+1
        //checking if we have a convergence
        if (checkConverged(i)) {
            copyFromSolution(solve,&E,&F,i);// copy it to solve
            //check if we have a convergence

            CalculateT(i,deltaT);//calculate Tn+1
            ApplyTandSource(i,deltaX,deltaT);
            return;
        }
        copyFromSolution(solve,&E,&F,i);// copy it to solve
        //check if we have a convergence

        CalculateT(i,deltaT);//calculate Tn+1
        ApplyTandSource(i,deltaX,deltaT);//we apply to solve Tn+1 and the src.
    }
}

void PredictorCorrectorSolution(int times,int i, void(*f)(),double deltaX,double deltaT) {

    PredictorCorrectorEF(1,i,f,deltaX,deltaT);
    //    CalculateT(i,deltaT);
    //ApplyTandSource(i,deltaX,deltaT);
    //after this step at E contains En+1 and T contains n+1 etc.
    //and not to forget solve contains En+1 and Fn+1
    //we apply to solve the addition of Tn+1 and the source

    //now we will want to solve the system again with Tn+1 instead of Tn
    //we copy the old E,F aka En and Fn to solve.
    /*    k = 0;
    int ok = 2;
    while (ok != 0) {
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
        ApplyTandSource(i,deltaX,deltaT);
        //and now we solve again
        PredictorCorrectorEF(2,i,f,deltaX,deltaT);
        //now we have the correct E and F. or the more correct. now we will re-calculate
        //T and apply the source again.
        CalculateT(i,deltaT);
        ApplyTandSource(i,deltaX,deltaT);
        ok--;
    }

    ApplyTandSource(i,deltaX,deltaT);*/
}

void CalculateT(int i,double deltaT) {
    int j;
    double sigma = 1.0;
    for ( j = 0; j < X; j++) {
      //this is where we calculate the next Temperture !
        T[j][i] = ((T[j][i-1] / deltaT) + sigma*E[j][i]) / (sigma + (1.0/deltaT) );
     }
}

void ApplyTandSource(int i,double deltaX,double deltaT) {
    int j;
    int k = 0;
    double Src;
    for ( j = 0; j < X; j++) {
        Src = 0.0;
        if (j*deltaX >= x0 || i*deltaT >= t0)
        {
            Src = 0;
        }
        //opacity here
         solve[2*j + 1] += getOpacity(j,timestep)*deltaT*T[j][i]+ Src*deltaT;
    }
    solve[0] = 0.5;
}

int checkConverged(int j) {
    int i;
    double resuSumE = 0.0;
    double resuSumF = 0.0;
    //solve contains En+1 this step, E[][] contains last step.
    for (i = 0; i < X; i++) {
        resuSumE += fabs(E[i][j] - solve[2*i + 1]);
    }
    for (i = 0; i <= X; i++) {
        resuSumF += fabs(F[i][j] - solve[2*i]);
    }
    if (resuSumE/X < 0.0000000001 && resuSumF/(X+1) < 0.0000000001) {
        return 1;
    }
    return 0;
}

double getOpacity(int space,int time1) {
        if (space >= X) {
          space = X-1;
        }
        double v = T[space][time1];
        double t  = pow(v,0.25);
        double a = 1.0/(pow(t,3.0) + 1e-15);
        return a;
}