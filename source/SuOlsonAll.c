#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "SuOlsonAll.h"
#include <malloc.h>
#include <string.h>
#include "tridFunc.h"
#define N 5000
#define X 5000
//#define NN (((X*2) + 1))
#define NN 5000
//#define N 10
//#define X 10
//#define NN (((X*2) + 1))
//#define NN 10
double epsilon = 0.00000000000000000000001;

void buildABLambdaT(int XX,int NX ,double (*EF)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void copyFromSolutionP1(double*solve,double(*mat)[X][N],double (*m)[X+1][N],int j);
void copyFromSolutionDiff(double*solve,double(*mat)[X][N],double (*m)[X+1][N],int j);
void sendToFileE(int);
void sendToFileT(int);
void PredictorCorrectorSolution(int times,int i,void(*f)(),void(*a)(),void(*)(),void(*)());
void constructLUDP1(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDDiff(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDP1AB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDP1MUAB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDDiffMUB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void CalculateT(int i,double deltaT);
void ApplyTandSourceP1(int i,double deltaX,double deltaT);
void ApplyTandSourceDiff(int i,double deltaX,double deltaT);
int checkConverged(int i);
double calculateMu(int space,int time1);
double calculateKappa(double);
double calculateA(double weff);
int setUpProgram(int,char*arg[]);
double getSource();
void setUpInitialCondition();

int currentTimeStep = 0;
double E[X][N],T[X][N];
//@@@CHANGED TO 2*NN +1
double L[NN],U[NN],mainD[NN],F[X+1][N],EF[X],D[X],solve[2*NN + 1];
double x0 = 0.5;
double t0 = 10.0;
double A[X];
double opacity;
double B[X];
double P1 = 0;
//double A = 3,B = 3;
double lambdaT;
double deltaX = 0.025;
double deltaT = 0.01;
int constOpacity = 0;
int main(int argc,char *argv[]) {
  double a,b,d;
  FILE*fp;
  int k,p,h,i=0,j=0;
  void (*funcptr) (int ,int ,double(*)[X],double[X][N],double[X+1][N],double[X][N],int);
  void (*BuildLUD)(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
  void (*applyTandS)(int,double,double);
  void (*copySolve)(double*,double(*mat)[X][N],double (*m)[X+1][N],int j);
  // now we have a trid-matrix size 2N on 2N
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR TEperture ie j const, i not. you move in X axis
  p = setUpProgram(argc,argv);
  //printf("0-Kershaw EF\n1-Levemore Pomraning EF\n2-Minerbo EF\n"
  //"3-P1\n4-P1AB\n5-Kershaw FL\n6-Levermore Pomraning FL\n7-Minerbo FL\n"
  //"8-Diffusion\n9-P1MUAB\n11-Asymptotic Diffusion\n12- Asymptotic Disc Diffusion\n");
  if (p != 100000) {
      if (!p) {
          funcptr = &buildEFKershaw;
          BuildLUD = &constructLUDP1;
          applyTandS = &ApplyTandSourceP1;
          copySolve = &copyFromSolutionP1;
      } else if (p == 1) {
          funcptr = &buildEFLP;
          BuildLUD = &constructLUDP1;
          applyTandS = &ApplyTandSourceP1;
          copySolve = &copyFromSolutionP1;
          BuildZ();
       }
       else if (p == 2) {
          funcptr = &buildEFMinerbo;
          BuildLUD = &constructLUDP1;
          applyTandS = &ApplyTandSourceP1;
          copySolve = &copyFromSolutionP1;
          BuildZ();
      } else if (p == 3) {
          deltaX = 0.01;
        funcptr = &buildNoEF;
        applyTandS = &ApplyTandSourceP1;
        copySolve = &copyFromSolutionP1;
        BuildLUD = &constructLUDP1;
       } else if (p == 4) {
        deltaX = 0.01;
        funcptr = &buildABLambdaT;
        applyTandS = &ApplyTandSourceP1;
        copySolve = &copyFromSolutionP1;
        BuildLUD = &constructLUDP1AB;
    } else if (p == 5) {
        deltaX = 0.01;
        funcptr = &buildDKershaw;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    }else if (p == 6) {
        deltaX = 0.01;
        funcptr = &buildDLP;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    }else if (p == 7) {
        deltaX = 0.01;
        funcptr = &buildDMinerbo;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    } else if (p == 8) {
        deltaX = 0.01;
        funcptr = &buildDDiff;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    }
    else if (p == 9) {
        deltaX = 0.01;
        funcptr = &buildABLambdaT;
        applyTandS = &ApplyTandSourceP1;
        copySolve = &copyFromSolutionP1;
        BuildLUD = &constructLUDP1MUAB;
    }
    else if (p == 11) {
        deltaX = 0.01;
        funcptr = &buildDDiffAsym;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
        BuildLUD = &constructLUDDiff;
    }
    else if (p == 12) {
        deltaX = 0.01;
        funcptr = &buildDDiscDiffAsym;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
        BuildLUD = &constructLUDDiffMUB;
    }
  }
  //boundry of the source
  setUpInitialCondition();
  //setting up the matrices
  for (currentTimeStep = 1; currentTimeStep < N; currentTimeStep++) {
      PredictorCorrectorSolution(1,currentTimeStep, funcptr,BuildLUD,applyTandS,copySolve);
    }

  for ( i = 0; i < N; i++) {
        for ( j = 0; j < X; j++) {
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000) {
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
                    printf("%f\t",pow(T[j][i],0.25));
             //printf("%f\t",E[j][i]);
           //           printf("%f\t",pow(E[j][i],0.5));
            }
        }
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000)
        printf("\n");
    }
    //printMatrix(E);
    if (constOpacity == 1) {
        sendToFileE(p);
    }
   else {
       sendToFileT(p);
   }
  return 0;
}


//copies from solve to the matrix
void copyFromSolutionP1(double*solve,double(*E)[X][N],double(*F)[X+1][N],int j) {
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
}

void copyFromSolutionDiff(double*solve,double(*E)[X][N],double(*F)[X+1][N],int j) {
    int i = 0;
    for ( i = 0; i < X; i++) {
      (*E)[i][j] = solve[i];
     
    }
  //  printf("\n");
}

void sendToFileE(int p) {
    int i = 0,j;
    FILE*fp;
    if (p == 0) {
        fp = fopen("../data/SuOlsonEddingtonFactorKershaw.txt","w");
    } else if (p == 1) {
        fp = fopen("../data/SuOlsonEddingtonFactorLP.txt","w");
    } else if (p == 2) {
        fp = fopen("../data/SuOlsonEddingtonFactorMinerbo.txt","w");
    } else if (p == 3 ) {
      fp = fopen("../data/SuOlsonP1Data.txt","w");
    } else if (p == 4) {
      fp = fopen("../data/SuOlsonP1ABData.txt","w");
  } else if (p == 5) {
      fp = fopen("../data/SuOlsonFluxLimitersDataKershaw.txt","w");
  }else if ( p == 6 ) {
      fp = fopen("../data/SuOlsonFluxLimitersDataLP.txt","w");
  }else if ( p == 7) {
      fp = fopen("../data/SuOlsonFluxLimitersDataMinerbo.txt","w");
  } else if ( p == 8) {
      fp = fopen("../data/SuOlsonData.txt","w");
  } else if ( p == 9) {
     fp = fopen("../data/SuOlsonP1MUABData.txt","w");
 }else if ( p == 11) {
    fp = fopen("../data/SuOlsonDiffusionAsymptoticData.txt","w");
  }else if ( p == 12) {
   fp = fopen("../data/SuOlsonDiffusionDiscAsymptoticData.txt","w");
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

void sendToFileT(int p) {
    int i = 0,j;
    FILE*fp;
    if (p == 0) {
        fp = fopen("../data/Temp/SuOlsonEddingtonFactorKershaw.txt","w");
    } else if (p == 1) {
        fp = fopen("../data/Temp/SuOlsonEddingtonFactorLP.txt","w");
    } else if (p == 2) {
        fp = fopen("../data/Temp/SuOlsonEddingtonFactorMinerbo.txt","w");
    } else if (p == 3 ) {
      fp = fopen("../data/Temp/SuOlsonP1Data.txt","w+");
    } else if (p == 4) {
      fp = fopen("../data/Temp/SuOlsonP1ABData.txt","w");
  } else if (p == 5) {
      fp = fopen("../data/Temp/SuOlsonFluxLimitersDataKershaw.txt","w");
  }else if ( p == 6 ) {
      fp = fopen("../data/Temp/SuOlsonFluxLimitersDataLP.txt","w");
  }else if ( p == 7) {
      fp = fopen("../data/Temp/SuOlsonFluxLimitersDataMinerbo.txt","w");
  } else if ( p == 8) {
      fp = fopen("../data/Temp/SuOlsonData.txt","w");
  } else if ( p == 9) {
     fp = fopen("../data/Temp/SuOlsonP1MUABData.txt","w");
 }else if ( p == 11) {
    fp = fopen("../data/Temp/SuOlsonDiffusionAsymptoticData.txt","w");
  }else if ( p == 12) {
   fp = fopen("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt","w");
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
            fprintf(fp,"%f ",pow(T[i][j],0.25));
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
}

void buildABLambdaT(int TTT,int nothing ,double (*EF)[X], double E1[X][N],double F1[X+1][N],double T1[X][N],int j) {
  int i,k;
  for (i = 0; i < X; i++) {
      double weff = calculateWeff(i,j);
      A[i] = calculateA(weff);
      B[i] = calculateB(weff);
      if (A[i] < epsilon) {
          A[i] = epsilon;
      }
      (*EF[i]) = 1.0/A[i];
  }
}

void PredictorCorrectorSolution(int times,int i, void(*f)(),void(*BuildLUD)(),void(*ApplyTS)(),void(*copySolve)()) {
    int j,k=0,p=0;
    //double copySolution[NN];
    /*for (j = 0; j < NN; j++) {
        if (j % 2 != 0) {
            copySolution[j] = E[k][i-1];
            k++;
        } else {
            copySolution[j] = F[p][i-1];
            p++;
        }
    }*/
    //we first do the basis, where we calculate E*,F*
    (*f)(X,N,&EF,E,F,T,i-1);//build EF or FL
    (*BuildLUD)(&L,&U,&mainD); // build LUD
    solveTriagonal(NN,&solve,L,U,mainD); // solve
    //now that we solved u(x,t+1), we will copy it to E.
    (*copySolve)(solve,&E,&F,i); //copy solution
    //note, when we solve the real E(n+1) we need copySolution
    CalculateT(i,deltaT);//we calculate Tn+1
    (*ApplyTS)(i,deltaX,deltaT);//we apply to solve Tn+1 and the src for the next step
   /*    times--;
    while (times) {
        times++;
        for ( j = 0; j < NN; j++) {
            solve[j] = copySolution[j];//copySolution contains E(n),F(n)
        }
        ApplyTandSourceP1(i,deltaX,deltaT);
        (*funcptr)(X,N,&EF,E,F,i);
        (*BuildLUD)(&L,&U,&mainD);
        solveTriagonalP1(&solve,L,U,mainD);//solve En+1 and Fn+1
        //checking if we have a convergence
        if (checkConverged(i)) {
            copyFromSolution(solve,&E,&F,i);// copy it to solve
            //check if we have a convergence

            CalculateT(i,deltaT);//calculate Tn+1
            ApplyTandSourceP1(i,deltaX,deltaT);
            return;
        }
        copyFromSolution(solve,&E,&F,i);// copy it to solve
        //check if we have a convergence

        CalculateT(i,deltaT);//calculate Tn+1
        ApplyTandSourceP1(i,deltaX,deltaT);//we apply to solve Tn+1 and the src.
    }*/
}

void constructLUDP1(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double opacity;
    for (i = 0; i < NN-1; i++) {
      if (i % 2 == 1 ) {//build E part
          (*U)[i] = deltaT/deltaX;
      } else {//build F part
          if (i != NN-1 && i != 0) {
            (*U)[i] = (deltaT*EF[j])/(deltaX);
          }
        j++;
      }
    }
    (*L)[0] = (*L)[1] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i % 2 == 1 ) {//E
             (*L)[i] = -deltaT/deltaX;
        } else {
            (*L)[i] = (-EF[j-1]*deltaT)/(deltaX);
            j++;
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
      if (i % 2 == 1) {//E
          opacity = getOpacity(i,currentTimeStep-1);
          (*mainD)[i] = 1.0 + deltaT*opacity;
      } else {//F
          opacity = getOpacity(i,currentTimeStep-1);
          (*mainD)[i] = 1.0 + deltaT*opacity;
          j++;
      }
    }
    //@@@@@@@added
    if (!constOpacity) {
      (*mainD)[0] = (double)1;
      (*U)[0] = -(double)1;
    }
}

void constructLUDP1AB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double opacity;
    for (i = 0; i < NN-1; i++) {
      if (i % 2 == 1) {//build E part
          (*U)[i] = deltaT/deltaX;
      } else {//build F part
          if (i != NN-1 && i != 0) {
            (*U)[i] = (deltaT)/(deltaX*A[j]);
          }
        j++;
      }
    }
    (*L)[0] = (*L)[1] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i % 2 == 1) {//E
             (*L)[i] = -deltaT/deltaX;
        } else {//F
            (*L)[i] = (-deltaT)/(deltaX*A[j]);
            j++;
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
      if (i % 2 == 1) {//E
          opacity = getOpacity(i,currentTimeStep-1);
          (*mainD)[i] = 1 + deltaT*opacity;
      } else {//F
           opacity = getOpacity(i,currentTimeStep-1);
          (*mainD)[i] = 1 + deltaT*opacity*B[j]/A[j];
          j++;
      }
    }
    //@@@ADDED
    if (!constOpacity) {
        (*mainD)[0] = (double)1.0;
        (*U)[0] = calculateMu(0,currentTimeStep-1);
    }
}

void constructLUDDiff(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double opacity;
    double lambda = deltaT/(deltaX*deltaX);
    for (i = 0; i < NN-1; i++) {
        if (i != X-1) {
            (*U)[i] = -lambda * (EF[i] + EF[i+1])/2.0;
        }
    }

    (*L)[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i != 0 ) {
            (*L)[i] = -lambda * (EF[i-1]+EF[i])/2.0;
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
        opacity = getOpacity(i,currentTimeStep-1);
        (*mainD)[i] = 1 + deltaT*opacity + lambda*(EF[i+1]+2*EF[i] + EF[i-1])/(2.0);
    }
    //@@@ added 2*deltaX
    double bb = 0;
    if (constOpacity == 0 ) {
        bb = deltaX;
    }
    (*mainD)[0] = 1 + deltaT*getOpacity(0,currentTimeStep-1) + lambda*(bb + EF[0] + EF[1])/2.0;
    i = NN - 1;
   // printf("%lf\t",EF[0] );
    opacity = getOpacity(i,currentTimeStep-1);
    (*mainD)[NN-1] = 1 + deltaT*opacity + lambda*(EF[i]+2*EF[i] + EF[i-1])/(2.0);
}

void constructLUDP1MUAB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    constructLUDP1AB(L,U,mainD);
    (*L)[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i % 2 == 1) {//E
             (*L)[i] = -deltaT/deltaX;
        } else {
            //(*L)[i] = (-EF[j-1]*deltaT)/(deltaX);
            double muprev,mucurrent,mu;
            muprev = calculateMu(j-1,currentTimeStep-1);
            mucurrent = calculateMu(j,currentTimeStep-1);
            mu = muprev/mucurrent;
            (*L)[i] = (-deltaT*mu)/(deltaX*A[j]);
            j++;
        }
    }
}

void constructLUDDiffMUB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double lambda = deltaT/(deltaX*deltaX);
    double mu,opacity;

    for (i = 0; i < NN-1; i++) {
        if (i != X-1) {
            mu = calculateMu(i+1,currentTimeStep-1);
            (*U)[i] = -lambda *mu*(EF[i] + EF[i+1])/2.0;
        }
    }

    (*L)[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i != 0 ) {
            mu = calculateMu(i-1,currentTimeStep-1);
            (*L)[i] = -lambda * mu*(EF[i-1]+EF[i])/2.0;
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
        mu = calculateMu(i,currentTimeStep-1);
        opacity = getOpacity(i,currentTimeStep-1);
        (*mainD)[i] = 1 + deltaT*opacity + lambda*mu*(EF[i+1]+2*EF[i] + EF[i-1])/(2.0);
    }
    //@@@ added 2*deltaX
    double bb= 0;
    if (constOpacity == 0 ) {
        bb = 2*deltaX*calculateMu(0,currentTimeStep-1);;
    }
    mu = calculateMu(0,currentTimeStep-1);
    opacity = getOpacity(0,currentTimeStep-1);
    (*mainD)[0] = 1 + deltaT*opacity + lambda*mu*(bb + EF[0] + EF[1])/2.0;
    i = NN - 1;
      mu = calculateMu(i,currentTimeStep-1);
      opacity = getOpacity(i,currentTimeStep-1);
      (*mainD)[NN-1] = 1 + deltaT*opacity + lambda*mu*(EF[i]+2*EF[i] + EF[i-1])/(2.0);
}

void CalculateT(int i,double deltaT) {
    int j;
    double sigmaA = 1.0;
    for ( j = 0; j < X; j++) {
        sigmaA = 1.0;
        // sigmaA = getOpacity(j,i-1);
         // printf("%lf\t",E[j][i]);
        T[j][i] = ((T[j][i-1] / deltaT) + sigmaA*E[j][i]) / (sigmaA + (1.0/deltaT) );
      //   printf("%lf\t",T[j][i]);
     }
    // printf("\n");
}

void ApplyTandSourceP1(int i,double deltaX,double deltaT) {
    int j;
    int k = 0;
    double Src;
    for ( j = 0; j < X; j++) {
        Src = 0;
        if (! (j*deltaX >= x0 || i*deltaT >= t0))
        {
            Src = getSource();
        }
         solve[2*j + 1] += getOpacity(j,i)*deltaT*T[j][i]+ Src*deltaT;
    }
    //@@@added this
    if (!constOpacity){
      solve[0] = (double)1/2;
      }

}

void ApplyTandSourceDiff(int i,double deltaX,double deltaT) {
    int j;
    int k = 0;
    double Src;
    for ( j = 0; j < X; j++) {
        Src = 0;
        if (! (j*deltaX >= x0 || i*deltaT >= t0))
        {
            Src = getSource();
        }
        solve[j] += getOpacity(j,i)*deltaT*T[j][i]+ Src*deltaT;
    }
    if (!constOpacity) {
        double l = deltaT / (deltaX*deltaX);
        solve[0] += l*deltaX/2.0;
      //  printf("%lf\n",solve[0]);
    }
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

double getdx() {
    return deltaX;
}

double getdt() {
    return deltaT;
}

double calculateWeff(int space,int time1) {
    double weff;
    double tt,ee;
    double opacity = getOpacity(space,time1);
    double Src = 0.0;
    if (! (time1*deltaT >= t0 || space*deltaX >= x0)) {
      Src = getSource();
    }
    tt = T[space][time1];
    ee = E[space][time1];
    weff = (opacity * tt + Src ) / (opacity*ee + 1e-15);
    return weff;
}

double calculateKappa(double weff) {
    if (weff < 0.01) {
        return 1;
    } else if (0.01 <= weff && weff <= 0.45) {
        double a = exp(-2.0/(weff + 1e-15));
        double b = (24.0 + 20.0*weff + 3.0*pow(weff,2.0)) / (pow(weff,2.0) + 1e-15);
        double c = (1.0- (4.0*a*(1.0 + a*((4.0 - 2.0*weff)/(weff + 1e-15)) + b*exp(-4.0/(weff + 1e-15)) ) ));
        return c;
    } else if (0.45 < weff && weff < 1.0) {
        return (1-weff) * calculateB(weff);
    } else  {
        return (weff-1)*calculateB(weff);
    }
}

double calculateB(double weff) {
    //double b1;
    if ( 0.59 <= weff && weff <=0.61) {
     return 1.0 / (0.80054 - 0.523*weff);
    } else {
        double c = (1.0 + weff) / 0.40528473;
        double a = (0.1326495 + weff*(0.03424169 + weff*(0.1774006 - weff)));
        double b = (0.3267567 + weff*(0.1587312 - weff*(0.5665676 + weff)));
        return (a*c)/b;
    }
}

double calculateA(double weff) {
    if (0.55 <= weff && weff <= 0.65) {
      return 0.96835 - 0.437*weff;
    } else {
        double a = ( 0.247 * (0.433 + 0.421*weff -2.681*weff*weff
             - 1.82*pow(weff,3.0) + 4.9*pow(weff,4.0) -1.06*pow(weff,5.0)
              + 2.56*pow(weff,6.0) ) );
        double b =  pow(0.33 + 0.159*weff - 0.567 * pow(weff,2.0) - pow(weff,3.0)  ,2.0);
      return a/b ;
    }
}

double calculateMu(int space,int time1) {
    double weff = calculateWeff(space,time1);
    double kappa = calculateKappa(weff);
    if (weff < 0.01) {
        return 1;
    } else if(0.01 <= weff && weff < 0.999) {
        if (kappa > 0) {
            double a = -weff/(2.0*kappa);
            double b = a*(log(1-kappa + 1e-15));
            return b;
        } else {
            return 1;
        }
    } else if (0.999<= weff && weff <= 1.001) {
        double a = 8.3548 + 1.5708 + weff;
        double b = 2.1228 + 2.4674*weff;
        return log(a/b);
    } else {
        double a = weff/(2*kappa + 1e-20);
        return a*log(1+kappa);
    }
}

double calculateMu2(double weff) {
    double kappa = calculateKappa(weff);
    if (weff < 0.01) {
        return 1;
    } else if(0.01 <= weff && weff < 0.999) {
        if (kappa > 0) {
            double a = -weff/(2.0*kappa);
            double b = a*(log(1-kappa + 1e-15));
            return b;
        } else {
            return 1;
        }
    } else if (0.999<= weff && weff <= 1.001) {
        double a = 8.3548 + 1.5708 + weff;
        double b = 2.1228 + 2.4674*weff;
        return log(a/b);
    } else {
        double a = weff/(2*kappa + 1e-20);
        return a*log(1+kappa);
    }
}

double getOpacity(int space,int time1) {
    if (constOpacity) {
        return 1;
    } else {
        double v = T[space][time1];
        double t  = pow(v,0.25);
        double a = 1.0/(pow(t,3.0) + 1e-15);;
        return a;
    }
}

int setUpProgram(int argc,char *argv[]) {
    //the first one is going to be if it is a constant opacity
    FILE*fp;
    int p = 0;
    char * line = NULL;
    fp = fopen(argv[1],"r");
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line,&len,fp)) != -1) {
    //    printf("%s\n",line);

        if (strstr(line,"Opacity:") != NULL) {
            int i = 1;
            while (line[i] != ':') {
                i++;
            }
            constOpacity = line[i+2] - '0';
        }

        else if (strstr(line,"Diffusion:") != NULL) {
            //this a diffusion type
            if (strstr(line,"Disc") != NULL) {
                p = 12;
                continue;
            } else if (strstr(line,"Asymptotic") != NULL) {
                p = 11;
                continue;
            } else if (strstr(line,"FL") != NULL) {
                if (strstr(line,"Kershaw") != NULL) {
                    p = 5;
                    continue;
                } else if (strstr(line,"LP") != NULL) {
                    p = 6;
                    continue;
                } else if (strstr(line,"Minerbo") != NULL) {
                    p = 7;
                    continue;
                }
            }
            p = 8;
        }

        else if (strstr(line,"P1: ") != NULL) {
            if (strstr(line,"MUAB") != NULL) {
                p = 9;
                P1 = 1;
                continue;
            } else if (strstr(line,"AB") != NULL) {
                p = 4;
                P1 = 1;
                continue;
            } else if (strstr(line,"EF") != NULL) {
                if (strstr(line,"Kershaw") != NULL) {
                    p = 0;P1 = 1;
                    continue;
                } else if (strstr(line,"LP") != NULL) {
                    p = 1;P1 = 1;
                    continue;
                } else if (strstr(line,"Minerbo") != NULL) {
                    p = 2;P1 = 1;
                    continue;
                }
            }
            p = 3;P1 = 1;
        }
    }
    printf("%d\n",p);
    fclose(fp);
    free(line);
    return p;
}

double getSource() {
    //we are in a const opacity i.e src is 1
    if (constOpacity) {
        return 1.0;
    } else {
        return 0;
    }
}

void setUpInitialCondition() {
    double Src;
    int i,j;
    for ( i = 0; i < NN; i++) {
        solve[i] = 0;
    }

    if (P1) {
      for ( i = 0; i < X; i++) {
        Src = 0;
        if (! (i*deltaX) >= x0) {
          Src = getSource();
        }
        solve[2*i + 1] = Src*deltaT;
      }
      if (!constOpacity) {
        solve[0] = (double)1/2;
      }
    }
    else {//diffusion
      for ( i = 0; i < X; i++) {
          Src = 0;
          if (! ((i)*deltaX >= x0)) {
              Src = getSource();
          }
          //@@@ add if & else
          if (i == 0 && constOpacity == 0) {
                  solve[i] = deltaT/(2.0*deltaX);
          } else {
              solve[i] = Src*deltaT;
          }
    }
  }
  //#pragma omp parallel for default(shared)
    if (constOpacity) {
        for (i = 0; i < X; i++) {
          for ( j = 0; j < N; j++) {
              T[i][j] = E[i][j] = 0.0;
          }
        }
    } else {
        for (i = 0; i < X; i++) {
          for ( j = 0; j < N; j++) {
              T[i][j] = E[i][j] = pow(10,-5);
          }
        }
    }
    for (i = 0; i < X+1; i++) {
      for ( j = 0; j < N; j++) {
          F[i][j] = 0.0;
      }
    }

    for ( i = 0; i < X; i++) {
        D[i] = EF[i] = (double)1.0/3.0;
    }
    for ( i = 0; i < X; i++) {
      A[i] = B[i] = 3.0;
    }

}
