
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "SuOlsonAll.h"
#include <malloc.h>
#include <string.h>
#include "nrutil.h"
#include "nrutil.c"
#include "tridag.c"
#include "tridFunc.h"
#define N 2000
#define X 2000
//#define NN (((X*2) + 1))
//#define NN 3001
//#define N 10
//#define X 10
#define NN 2000
//#define NN X
//#define NN 10
double epsilon = 1e-20;
void buildABLambdaT(int XX,int NX ,double (*EF)[X], double E[X + 1][N],double F[X][N],double [X][N],int j);
void copyFromSolutionP1(double*solve,double(*mat)[X + 1][N],double (*m)[X][N],int j);
void copyFromSolutionDiff(double*solve,double(*mat)[X + 1][N],double (*m)[X][N],int j);
void sendToFileE(int);
void sendToFileT(int);
void PredictorCorrectorSolution(int times,int i,void(*f)(),void(*a)(),void(*)(),void(*)());
void constructLUDP1(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDDiff(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDP1AB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDP1MUAB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void constructLUDDiffMUB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
void sendToFileW(int p);
void CalculateT(int i,double deltaT);
void calculateFlux(int i);
void ApplyTandSourceP1(int i,double deltaX,double deltaT);
void ApplyTandSourceDiff(int i,double deltaX,double deltaT);
void findWavefront(int);
inline double getInitValue();
int checkConverged(int i);
double calculateMu(int space,int time1);
double calculateKappa(double);
double calculateA(double weff);
int setUpProgram(int,char*arg[]);
double getSource(int space,int time1);
double getT(int space,int time1);
double getT2(int space, int time1);
void setUpInitialCondition();
double getTH();
double Avg1(double,double);
double convertLineToDouble(char* str, int len);
int convertLineToInt(char* str, int len);
double Avg2(double,double);
void update_dt();
void applyBC(int);
double Min(double, double);
int currentTimeStep = 0;
double F[X][N],T[X + 1][N];
//@@@CHANGED TO 2*NN +1
double L[NN + 1],U[NN + 1],mainD[NN + 1],E[X][N],EF[X],D[X],solve[2 * NN + 1],Weff[X+1][N];
double waveFront[X];
double currentTime[NN];
int problem;
double initV = 0.000000000000001;
double TH[2][759];
double a_[NN + 1], b_[NN + 1], si_[NN + 1], r_[NN + 1];
double x0 = 0.5;
double t0 = 10.0;
double A[X+1];
double eps = 1.0;
double opacity;
double B[X+1];
double sigma_boltzman = 5.670373e-5;
double c = 3E10;
// 50 mg/cm^3
//double rhoSilicon = 50;
//convert it to g/cm^3
double rho = 0.05;
//double rho = 1.0;
//double c = 1;
double arad = 7.56E-15;
double Cv = 0;
double d_frac = 0.1;
double P1 = 0;
//-silicon is 3.53
double alpha = 3.5;
double beta = 1.1;
//double alpha = 3.0;
//double beta = 1.0;
double mu_sio = 0.1;
//double mu_sio = 1.0;
double s_lambda = 0.75;
double s_f = 8.8E13;
//double s_f =  4.0 * 7.56E-15;
//in g/cm^3
double s_g = 1.0/9175.0;
//double s_g = 1.0;
//double alpha = 1.0;
int Classic = 0;
//double A = 3,B = 3;
double lambdaT;
double bb = 0;
double previousWeff = 0;
//in cm
//double deltaX = 0.01;
//in mm it is 0.001
//in seconds
//double deltaT = 0.01;
double deltaT = 0.01;
double deltaX = 0.01;
double th;
int constOpacity = 0;

int main(int argc,char *argv[]) {
  double a,b,d;
  FILE*fp;
  int k,p,h,i=0,j=0;
  void (*funcptr) (int ,int ,double(*)[X],double[X][N],double[X+1][N],double[X][N],int);
  void (*BuildLUD)(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]);
  void (*applyTandS)(int,double,double);
  void (*copySolve)(double*,double(*mat)[X+1][N],double (*m)[X][N],int j);
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
        funcptr = &buildNoEF;
        applyTandS = &ApplyTandSourceP1;
        copySolve = &copyFromSolutionP1;
        BuildLUD = &constructLUDP1;
       } else if (p == 4) {
        funcptr = &buildABLambdaT;
        applyTandS = &ApplyTandSourceP1;
        copySolve = &copyFromSolutionP1;
        BuildLUD = &constructLUDP1AB;
    } else if (p == 5) {
        funcptr = &buildDKershaw;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    }else if (p == 6) {
        funcptr = &buildDLP;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    }else if (p == 7) {
        funcptr = &buildDMinerbo;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    } else if (p == 8) {
        funcptr = &buildDDiff;
        BuildLUD = &constructLUDDiff;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
    }
    else if (p == 9) {
        //deltaX = 0.01;
        funcptr = &buildABLambdaT;
        applyTandS = &ApplyTandSourceP1;
        copySolve = &copyFromSolutionP1;
        BuildLUD = &constructLUDP1MUAB;
    }
    else if (p == 11) {
        funcptr = &buildDDiffAsym;
        applyTandS = &ApplyTandSourceDiff;
        copySolve = &copyFromSolutionDiff;
        BuildLUD = &constructLUDDiff;
    }
    else if (p == 12) {
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
    printf("%d\n",constOpacity);
  for ( i = 0; i < N; i++) {
        for ( j = 0; j < X; j++) {
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000) {
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
               //     printf("%f\t",getT(j,i));
            printf("%f\t",E[j][i]);
           //           printf("%f\t",pow(E[j][i],0.5));
            }
        }
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000)
        printf("\n");
    }
   // printf("\n");
    if (constOpacity == 1) {
        sendToFileE(p);
    }
   else {
        sendToFileE(p);
        sendToFileT(p);
        sendToFileW(p);
   }
  return 0;
}

//copies from solve to the matrix
void copyFromSolutionP1(double*solve,double(*E)[X + 1][N],double(*F)[X][N],int j) {
  int i = 0,k = 0;
  for ( i = 0; i < NN; i+=2) {
      if (k == X) {
          return;
      }
    (*F)[k][j] = solve[i+1];
    (*E)[k][j] = solve[i];
    k++;
  }
  (*E)[X][j] = solve[NN-1];
}

void copyFromSolutionDiff(double*solve,double(*E)[X + 1][N],double(*F)[X][N],int j) {
    int i = 0;
    for ( i = 0; i < X; i++) {
      (*E)[i][j] = solve[i];
    }
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
      fprintf(fp, "%f ", currentTime[i] * c );
    }
    fprintf(fp, "\n");
    for ( j = 0; j < N; j++) {
      for ( i = 0; i < X; i++) {
            fprintf(fp,"%f ", (E[i][j]));
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
 } else if ( p == 11) {
    fp = fopen("../data/Temp/SuOlsonDiffusionAsymptoticData.txt","w");
  } else if ( p == 12) {
   fp = fopen("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt","w");
  }
    for ( i = 0; i < N; i++) {
            fprintf(fp, "%10e ",deltaX*(i) );
    }
    fprintf(fp, "\n");
    double kk = c;
    if (problem == 2){
        kk = 1E9;
    }
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%10e ",currentTime[i] * kk );
    }
    fprintf(fp, "\n");
    if (problem != 2){
        for ( j = 0; j < N; j++) {
            for ( i = 0; i < N; i++) {
                fprintf(fp,"%10e ",(getT(i,j)));
            }
        fprintf(fp,"\n");
        }
    } else {
        for ( j = 0; j < N; j++) {
            for ( i = 0; i < N; i++) {
                fprintf(fp,"%f ",(getT(i,j)/11605.0));
            }
        fprintf(fp,"\n");
        }
    }
    fclose(fp);
}

void sendToFileW(int p) {
    int i = 0,j;
    FILE*fp;
   if ( p == 8) {
      fp = fopen("../data/Temp/Back_1500_WaveFront.txt","w");
   }
    for ( i = 0; i < N; i++) {
        fprintf(fp, "%10e ",deltaX*(i) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%10e ",currentTime[i] * 1E9 );
    }
    fprintf(fp, "\n");
    for ( j = 0; j < X; j++) {
        fprintf(fp,"%10e ",waveFront[j]);
    
    }
      fprintf(fp,"\n");
      for ( j = 0; j < X; j++) {
        fprintf(fp,"%10e ",sigma_boltzman * E[0][j]/arad);
    
    }
      fprintf(fp,"\n");
    for ( j = 0; j < X; j++) {
        fprintf(fp,"%f ", F[0][j]/2.0);
    
    }
      fprintf(fp,"\n");
    fclose(fp);
}

void buildABLambdaT(int TTT,int nothing ,double (*EF)[X], double E1[X + 1][N],double F1[X][N],double T1[X][N],int j) {
  int i,k;
  return;
  for (i = 0; i < X + 1; i++) {
      double weff = calculateWeff(i,j);
      A[i] = calculateA(weff);
      B[i] = calculateB(weff);
      if (A[i] < epsilon) {
          A[i] = epsilon;
      }
  }
  A[X] = A[X-1];
  B[X] = B[X-1];
}

void PredictorCorrectorSolution(int times,int i, void(*f)(),void(*BuildLUD)(),void(*ApplyTS)(),void(*copySolve)()) {
    int j,k=0,p=0;
    double solve_prev[N + 1];
    //we first do the basis, where we calculate E*,F*
    for (j = 0; j < X; j++) {
        solve_prev[j] = solve[j];
    }
    ApplyTandSourceDiff(i - 1, deltaX, deltaT);
    applyBC(i - 1);
    (*f)(X,N,&EF,E,F,T,i - 1);//build EF or FL
    (*BuildLUD)(&L,&U,&mainD); // build LUD
    solveTriagonal(NN, &solve, L, U, mainD);
    (*copySolve)(solve,&E,&F,i);
    CalculateT(i, deltaT);
    update_dt();
    return;
    //for second predictor corrector..
    for (j = 0; j < X; j++) {
        solve[j] = solve_prev[j];
    }
    ApplyTandSourceDiff(i, deltaX, deltaT);
    applyBC(i);
    (*f)(X,N,&EF,E,F,T,i);//build EF or FL
    (*BuildLUD)(&L,&U,&mainD); // build LUD
    solveTriagonal(NN, &solve, L, U, mainD);
    (*copySolve)(solve,&E,&F,i);
    CalculateT(i, deltaT);
   // findWavefront(i - 1);
    //calculateFlux(i - 1);
    
    return;

}

void constructLUDP1(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double opacity;
    for (i = 0; i < NN-1; i++) {
      if (i % 2 == 0) {//build E part
            (*U)[i] = deltaT/(deltaX);
      } else {//build F part
            (*U)[i] = c*c*deltaT/(deltaX);
      }
    }

    (*L)[0] = (*L)[1] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i % 2 == 0 ) {//E
               (*L)[i] = -1.0*deltaT/(deltaX);
        } else {
               (*L)[i] = (-c*c*deltaT)/(deltaX);
        }
    }
    j = 0;
    int k = 0;
    for (i = 0; i < NN; i++) {
      if (i % 2 == 0) {//E
          opacity = getOpacity(k,currentTimeStep-1);
          (*mainD)[i] = 1.0 + opacity*deltaT*c;
          k++;
      } else {//F
          double opacityprev = getOpacity(k-1,currentTimeStep-1);
          double opacitycurr = getOpacity(k,currentTimeStep-1);
          opacity = getOpacity(k - 1,currentTimeStep-1);
          opacity = Avg2(opacitycurr,opacityprev);
          (*mainD)[i] = 3.0 + 3.0*opacity*deltaT*c;
          j++;
      }
    }
}

void constructLUDP1AB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double opacity;
    double bAvg,muAvg,aAvg;
    for (i = 0; i < NN-1; i++) {
      if (i % 2 == 0 ) {//build E part
            (*U)[i] = deltaT/(deltaX);
      } else {//build F part
            (*U)[i] = (c*c*deltaT)/(deltaX);
      }
    }
    
    (*L)[0] = (*L)[1] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i % 2 == 0) {//E
            (*L)[i] = -deltaT/(deltaX);
        } else {//F
            (*L)[i] = (-c*deltaT*c)/(deltaX);
        }
    }
    j = 0;
    int k = 0;
    for (i = 0; i < NN; i++) {
      if (i % 2 == 0) {//E
          opacity = getOpacity(k,currentTimeStep-1);
          (*mainD)[i] = 1.0 + opacity*deltaT*c;
           k++;
      } else {//F
          double opacityprev = getOpacity(k-1,currentTimeStep-1);
          double opacitycurr = getOpacity(k,currentTimeStep-1);
          opacity = Avg2(opacitycurr,opacityprev);
          (*mainD)[i] = A[j] + opacity*B[j]*deltaT*c;
          j++;
      }
    }

}

void constructLUDP1MUAB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    constructLUDP1AB(L,U,mainD);
     //(*L)[0] = 0.0;
    j = 0;
    int k = 0;
    for ( i = 0; i < NN; i++) {
        if (i % 2 == 0) {//E
            k++;
        } else {
            double muprev,mucurrent;
            double bAvg,muAvg,aAvg;
            muprev = calculateMu2(calculateWeffNonAvg(j, 
            currentTimeStep - 1));
            mucurrent = calculateMu2(calculateWeffNonAvg(j , 
            currentTimeStep - 1));
            double udiv = mucurrent;
            double ldiv1 =  muprev;
            
            (*U)[i] = (c*udiv*c*deltaT)/(deltaX);
            (*L)[i] = -(c*ldiv1*c*deltaT)/(deltaX);
            double opacityprev = getOpacity(k - 1,currentTimeStep-1);
            double opacitycurr = getOpacity(k,currentTimeStep-1);
            opacity = Avg2(opacitycurr,opacityprev);
            double mu = calculateMu(j , currentTimeStep - 1);
            (*mainD)[i] = A[j]*mu + mu*opacity*B[j]*(c*deltaT);
            if (currentTimeStep == 273) {
               // printf("%lf\t%lf\t%lf\t%d\n",udiv,ldiv1,mu,i/2);
               
            }
            j++;
        }
    }
    (*L)[0] = 0.0;

}

void constructLUDDiff(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
   int i,j = 0;
    double opacity;

    double lambda = deltaT/(deltaX*deltaX);
    for (i = 0; i < NN-1; i++) {
        if (i != X-1) {
            (*U)[i] = -lambda *c* (EF[i] + EF[i+1])/2.0;
        }
    }

    (*L)[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i != 0 ) {
            (*L)[i] = -lambda * c * (EF[i -1] + EF[i])/2.0;
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
        opacity = getOpacity(i, currentTimeStep - 1);
        (*mainD)[i] = 1 + deltaT*c*opacity + lambda*c*(EF[1 + i] + 2*EF[i] + EF[i - 1])/(2.0);
    }
    //@@@ added 2*deltaX
   // double bb = 0;
    if (constOpacity != 1 ) {
       // bb = -2*pow(EF[0], 2) / (EF[0] + 0.5*deltaX); 
    }
    (*mainD)[0] = 1.0 + deltaT* c * getOpacity(0, currentTimeStep-1) + lambda*c*(bb + EF[0] + EF[1])/2.0;
    i = NN - 1;
    if (constOpacity == -1){
       // (*mainD)[0] = 1.0;
       // (*U)[0] = 0.0;
    }
   // printf("%lf\t",EF[0] );
    opacity = getOpacity(i,currentTimeStep-1);
    (*mainD)[NN-1] = 1.0 + deltaT*c*opacity + lambda* c *(EF[i]+2*EF[i] + EF[i-1])/(2.0);
}

void constructLUDDiffMUB(double (*L)[NN],double (*U)[NN],double (*mainD)[NN]) {
    int i,j = 0;
    double lambda = deltaT/(deltaX*deltaX);
    double mu,opacity;

    for (i = 0; i < NN-1; i++) {
        if (i != X-1) {
            mu = calculateMu(i+1,currentTimeStep-1);
          //      mu = calculateMu2(calculateWeffNonAvg(i+1,currentTimeStep-1));
            (*U)[i] = -lambda * c *mu*(EF[i] + EF[i+1])/2.0;
        }
    }

    (*L)[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i != 0 ) {
            mu = calculateMu(i-1,currentTimeStep-1);
          //   mu = calculateMu2(calculateWeffNonAvg(i-1,currentTimeStep-1));
            (*L)[i] = -lambda *c* mu*(EF[i-1]+EF[i])/2.0;
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
        mu = calculateMu(i,currentTimeStep-1);
      //  mu = calculateMu2(calculateWeffNonAvg(i,currentTimeStep-1));
        opacity = getOpacity(i,currentTimeStep-1);
      (*mainD)[i] = 1 + deltaT*c*opacity + lambda*c*mu*(EF[i+1]+2*EF[i] + EF[i-1])/(2.0);
    }
    //@@@ added 2*deltaX
    double bb= 0;
    if (constOpacity == 0 ) {
        bb = 2.0*deltaX*calculateMu(0,currentTimeStep-1);
    }
    mu = calculateMu(0,currentTimeStep-1);
    opacity = getOpacity(0,currentTimeStep-1);
    (*mainD)[0] = 1.0 + deltaT*c*opacity + lambda*c*mu*(bb + EF[0] + EF[1])/2.0;
    i = NN - 1;
      mu = calculateMu(i,currentTimeStep-1);
      opacity = getOpacity(i,currentTimeStep-1);
      (*mainD)[NN-1] = 1 + deltaT*c*opacity + lambda*c*mu*(EF[i]+2*EF[i] + EF[i-1])/(2.0);
}

void CalculateT(int i,double deltaT) {
    int j;
    int stop = 3;
    for ( j = 0; j < X; j++) {
        double t = getT(j, i - 1);
        double cap = getCv(j, i - 1);
        double coeff = (getOpacity(j, i - 1) * 4.0  * pow(t, 3) * arad)
        / (cap);
        //if (coeff != 1.0)
           // coeff = 1.0;
           coeff = deltaT*c*coeff;
        T[j][i] = ((T[j][i - 1]/coeff) + (E[j][i])) /
                 (1.0 + 1.0/coeff);
        //T[j][i] = ((T[j][i-1] + deltaT * c * E[j][i]))/(1.0 + c*deltaT) ;

    }  
}

void ApplyTandSourceP1(int i,double deltaX,double deltaT) {
    int j;
    int k = 0;
    double Src;
    for ( j = 0; j < X  ; j++) {
        Src = getSource(j,i);
        printf("a");
        solve[2*j] += getOpacity(j,i)*deltaT*c*T[j][i]+ Src*deltaT*c;
    }
    //@@@added this
    if (!constOpacity) {
      applyBC(1);
    }
}

void ApplyTandSourceDiff(int i,double deltaX,double deltaT) {
    int j;
    int k = 0;
    double Src;
    for ( j = 0; j < X; j++) {
        Src = getSource(j, i);
        solve[j] += getOpacity(j, i)*deltaT*c*T[j][i] + Src*c*deltaT;
    }
    //applyBC(currentTimeStep - 1);
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
    double Src = getSource(space,time1);
    if (space >= X){
        space = X-1;
    }
    tt = T[space][time1];
    ee = E[space][time1];
    //weff for non-classic, avg
    if (!constOpacity && P1) {
        //avg tempreture
        double t1 = (getT(space, time1));
        double t2 = getT(space - 1,time1);
        if (t1 != t1 || t2 != t2){
        //    printf("6. Time:%d\tT1: %lf\t T2:%lf\n",space,t1,t2);
        }
        tt = pow(Avg1(t1, t2),4);
        tt = arad*tt;
        //avg energy
        ee = Avg2(E[space - 1][time1], E[space][time1]);
        //avg opacity
        opacity = getOpacity(space ,time1);
        double op1 = getOpacity(space-1 ,time1);
        opacity = Avg2(opacity, op1);
    }   
    if ( !constOpacity && space == 0) {
        ee = E[space][time1];
        tt = T[space][time1];
        opacity = getOpacity(space,time1);
            Src = (2.0*getFinc() - (c*ee/2.0)) /c ;
        }

    weff = (opacity * tt + Src ) / (opacity*ee + 1e-15);
    if (weff != weff || ee != ee  || opacity != opacity) {
        printf("2. Time:%d \tTempreture:%lf\tEnergy:%15.15lf\tWeff: %lf\n",currentTimeStep, tt,ee,weff);
        exit(0);
    }
    return weff;
}

double calculateWeffNonAvg(int space,int time1) {
    double weff;
    double tt,ee;
    double opacity = getOpacity(space,time1);
    double Src = getSource(space,time1);
    tt = T[space][time1];
    ee = E[space][time1];
    if (space >= X){
        space = X-1;
    }
    //weff for non-classic, avg
    if (!constOpacity) {
        if (space == 0) {
            Src = (2.0*getFinc() - c*E[0][time1]/2.0)/c;
        }
    }
    weff = (opacity * tt + Src ) / (opacity * ee /*+ 1e-15*/);
    if (weff != weff || ee != ee) {
        printf("1. Time:%d \tTempreture:%lf\n",currentTimeStep, tt);
        exit(0);
    }
    return weff;
}

double calculateKappa(double weff) {
    if (weff < 0.01) {
        return 1.0;
    } else if (0.01 <= weff && weff <= 0.45) {
        double aa = exp(-2.0/(weff));
        double bb = (24.0 + 20.0*weff + 3.0*pow(weff,2.0)) / (pow(weff,2.0) );
        double xyz = (1.0- (4.0*aa*(1.0 + aa*((4.0 - 2.0*weff)/(weff)) + bb*exp(-4.0/(weff)) ) ));
        return xyz;
    } else if (0.45 < weff && weff < 1.0) {
        return (1.0-weff) * calculateB(weff);
    } else  {
        return (weff-1.0) * calculateB(weff);
    }
}

double calculateB(double weff) {
    //double b1;
    
    if ( 0.59 <= weff && weff <=0.61) {
     return 1.0 / (0.80054 - 0.523*weff);
    } else {
        double xyz = (1.0 + weff) / 0.40528473;
        double aa = (0.1326495 + weff*(0.03424169 + weff*(0.1774006 - weff)));
        double ba = (0.3267567 + weff*(0.1587312 - weff*(0.5665676 + weff)));
        return (aa*xyz)/ba;
    }
}

double calculateA(double weff) {
    if (0.55 <= weff && weff <= 0.65) {
      return 0.96835 - 0.437*weff;
    } else {
        double aa = ( 0.247 * (0.433 + 0.421*weff -2.681*weff*weff
             - 1.82*pow(weff,3.0) + 4.9*pow(weff,4.0) -1.058*pow(weff,5.0)
              + 2.56*pow(weff,6.0) ) );
        double ba =  pow(0.327 + 0.159*weff - 0.567 * pow(weff,2.0) - pow(weff,3.0)  ,2.0);
      return aa/ba ;
    }
}

double calculateMu(int space,int time1) {
    double weff = calculateWeff(space,time1);
    return calculateMu2(weff);
    double kappa = calculateKappa(weff);
    if (weff < 0.01) {
        return 1.0;
    } else if(0.01 <= weff && weff < 0.999) {
        if (kappa > 0) {
            double a = -weff/(2.0*kappa);
            double b = a*(log(1.0-kappa + 1e-15));
            return b;
        } else {
            return 1;
        }
    } else if (0.999<= weff && weff <= 1.001) {
        double a = 8.3548 + 1.5708 + weff;
        double b = 2.1228 + 2.4674*weff;
        return log(a/b);
    } else {
        double a = weff/(2.0*kappa + 1e-15);
        return a*log(1.0+kappa);
    }
}

double calculateMu2(double weff) {
    double kappa = calculateKappa(weff);
    if (weff < 0.01) {
        return 1;
    } else if(0.01 <= weff && weff < 0.999) {
        if (kappa > 0) {
            double a = -weff/(2.0*kappa);
            double b = a*(log(1.0 - kappa + 1e-15));
            return b;
        } else {
            return 1.0;
        }
    } else if (0.999<= weff && weff <= 1.001) {
        double a = 8.3548 + 1.5708 + weff;
        double b = 2.1228 + 2.4674*weff;
        return log(a/b);
    } else {
        double a = weff/(2.0*kappa /*+ 1e-15*/);
        return a*log(1.0 + kappa);
    }
}

double getOpacity(int space, int time1) {
    if (constOpacity == 1) {
        return 1.0;
    } else if (constOpacity == 0){
        if (space >= X) {
          space = X-1;
        }
        double t = getT(space, time1);
        double a = 1.0/(pow(t, 3.0));
        return a;
    } else {
        double t_galpha = (pow(rho, s_lambda + 1) 
                          / (s_g * pow(getT(space, time1), alpha)));
        //double bbbb = rho / t_galpha; 
        double t = getT(space, time1);
        double a = 1.0/(pow(t, 3.0));
        if (t_galpha != a) {
           // printf("bad opacity\n");
            //exit(1);
        }
        return t_galpha;
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
        p = 8;
        // set up the constants of the problem
        if (strstr(line, "SuOlson") != NULL) {
            problem = 0;
        } else if (strstr(line, "Olson") != NULL) {
            problem = 1;
        } else if (strstr(line, "Back") != NULL) {
            problem = 2;
        } else if (strstr(line, "x0") != NULL) {
            x0 = convertLineToDouble(line, len);
        } else if(strstr(line, "t0") != NULL) {
            t0 = convertLineToDouble(line, len);
        } else if(strstr(line, "g") != NULL) {
            s_g = convertLineToDouble(line, len);
        } else if(strstr(line, "f") != NULL) {
            s_f = convertLineToDouble(line, len);
        } else if(strstr(line, "mu") != NULL) {
            mu_sio = convertLineToDouble(line, len);
        } else if(strstr(line, "lambda") != NULL) {
            s_lambda = convertLineToDouble(line, len);
        } else if(strstr(line, "deltaT") != NULL) {
            deltaT = convertLineToDouble(line, len);
        } else if(strstr(line, "deltaX") != NULL) {
            deltaX = convertLineToDouble(line, len);
        } else if(strstr(line, "alpha") != NULL) {
            alpha = convertLineToDouble(line, len);
        } else if(strstr(line, "beta") != NULL) {
            beta = convertLineToDouble(line, len);
        } else if(strstr(line, "rho") != NULL) {
            rho = convertLineToDouble(line, len);
        } else if(strstr(line, "TH") != NULL) {
            th = convertLineToDouble(line, len);
        } else if(strstr(line, "Temp_0") != NULL) {
            initV = convertLineToDouble(line, len);
        } else if(strstr(line, "dfrac") != NULL) {
            d_frac = convertLineToDouble(line, len);
        }

        else if (strstr(line,"Opacity:") != NULL) {
            int i = 1;
            while (line[i] != ':') {
                i++;
            }
            constOpacity = line[i+2] - '0';
            constOpacity = -1;
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
    if (problem == 0) {
        deltaT = deltaT / c; //units
        bb = 0;
    } else if (problem == 2) {
      //  deltaT = deltaT/c;
       deltaT = deltaT * 1E-9; //units
        bb = deltaX;
        s_f = s_f / (pow(1160452.0, beta));
        s_g = s_g / pow(1160452.0, alpha);
    } else if (problem == 1) {
        deltaT = deltaT / c;
        s_f = 4.0 * arad;
        bb = deltaX;
    }
    printf("%d\n",p);
    fclose(fp);
    free(line);
    return p;
}

double getSource(int space,int time1) {
    //we are in a const opacity i.e src is 1
         if ( (space*deltaX < x0 && currentTime[time1]*c < t0))
        {
            return 1.0;
        }
        return 0.0;
}

void setUpInitialCondition() {
    double Src;
    int i,j;
    //arad = 4.0 * sigma_boltzman / c;
    //deltaT = deltaT / c;
    //deltaT = deltaT / c;
    //s_f = s_f / (pow(1160452.0, beta));
    //s_g = s_g / pow(1160452.0, alpha);
    //s_f = 4.0 * arad;
    //s_g = 1.0;
    //rho = 1.0;
    //beta = 1.0;
    //alpha = 3.0;
   // alpha = 4.0 * arad;


    for ( j = 0; j < N; j++) {
        E[j][0] = T[j][0] = pow(initV, 4) * arad;
    }
   

    for ( i = 0; i < X; i++) {
        D[i] = EF[i] = (double)1.0/(3.0 * getOpacity(i, 0));
        solve[i] = E[i][0];
      //  solve[i] = getOpacity(i, 0) * deltaT * c * T[i][0] + E[i][0] +  getSource(i, 0) *c*deltaT;
       // printf("%10e\n",getOpacity(i, 0));
    }
    //solve[0] += deltaT*c*arad/(2.0*deltaX);
   // applyBC(0);
    for ( i = 0; i < X; i++) {
      A[i] = B[i] = 3.0;
    }
    FILE *fp1;
    char buff[255];
    
    fp1 = fopen("dataset3.csv","r");
    for ( i = 0; i < 190; i++) {
        fscanf(fp1,"%s",buff);
        TH[0][i] = atof(buff);
        fscanf(fp1,"%s",buff);
        TH[1][i] = atof(buff);
        if (TH[1][i] < 0 ){
            printf("negative");
            exit(1);
        }
    }
    currentTime[0] = deltaT;
    fclose(fp1);
}

double getCv(int space, int time1) {
    if (constOpacity == -1) {
       double abb = beta * s_f * 
        pow(getT(space, time1), beta - 1.0) * pow(rho,-mu_sio);
        //if (abb != 4.0*arad)
          //  printf("bad cv\t%10e\t%10e\n",s_f, 4.0*arad);
        return abb;
    } else if (constOpacity == 1) {
        double abb = getT(space,time1);
        return (alpha*pow(abb,3));
    } else {
        return 4.0 * arad * pow(getTH(time1), 3);
    }
}

double getT(int space,int time1) {
    return pow(T[space][time1] / arad,0.25);
}

double getT2(int space,int time1) {
    return pow(T[space][time1] /( arad),0.25);
}

double getTH(int time1){
    return th;
    int i,j;
    return 190.0 * 11604.52;
    double t = time1*deltaT*1E9;
    for (i = 1; i < 190; i++){
        //printf("%10e\t%10e\n",t,TH[0][i]);
        if (TH[0][i - 1] < t && t <  TH[0][i]){
            return 11604.52* TH[1][i];
        }
    }
}

double getFinc(){
    return (c*arad*pow(getTH(1) ,4)/4.0);
}

void applyBC(int time1) {
  if (constOpacity == 0) {
    //double mu = calculateMu2(calculateWeff(0,currentTimeStep ));
    //solve[0] += getFinc() * deltaT * 2.0/(c * deltaX);
   double l = deltaT / (deltaX*deltaX);
    solve[0] += l*c*arad*deltaX/2.0;
    //solve[0] += (2.0*getFinc() - c*solve[0]*mu)*(deltaT/deltaX);
  } else if (constOpacity == -1) {
      //maybe currenttimestep - 1...
     // printf("%15.15lf\t",solve[0]);
    solve[0] += arad * c * pow(getTH(time1) ,4) * deltaT/(2.0*deltaX);
  // solve[0] += (2.0*getFinc() * deltaT / deltaX*(EF[0] + 0.5*deltaX));
  //  solve[0] += ( (2.0*getFinc() *EF[0]* (deltaT/deltaX)) / (EF[0] + deltaX*0.5)); 
    //solve[0] = arad * pow(getTH(time1), 4);
    //  double l = deltaT / (deltaX*deltaX);
    //solve[0] += l*c*arad*deltaX/2.0;  
    //   solve[0] += getFinc() * deltaT * 2.0/(c * deltaX);
    //  printf("%15.15lf\n",solve[0]);
  } else if (constOpacity == 1) {
      // solve[0] += getFinc() * deltaT * 2.0/(c * deltaX);
  }
}

double Avg1(double xx,double yy) {
    return (xx + yy) / 2.0;
}

double Avg2(double xx,double yy) {
    if ( xx == 0.0 || yy == 0.0) {
        return Avg1(xx,yy);
    }
   return (2.0 * xx * yy)/(xx + yy);
}

inline double getInitValue() {
    if (constOpacity == -1) {
        return arad * pow(190 * 11604.505, 4);
    } else {
      double tm = pow(10,-1.25)*5*getTH(1);
      return arad*pow(tm,4);
      }
}

void findWavefront(int time1) {
    int i, j;
    double maxT = 0;
    double position;
    double tol = 1.0 - 0.1;
    for(i = 1; i < X; i++) {
        if (getT(i, time1) > maxT ) {
            maxT = getT(i, time1);
            position = i * deltaX;
          //  printf("%15.15lf\n",getT(i,time1));
        }
    }
  //  printf("%15.15lf\t%lf\t",maxT,position);
    tol = maxT * tol; // tolerance T value of 5% 
    //tol is the tolerance value..
    //we find the maximum value of T, but with a tolerance of 5% from the real maximum value.
    maxT = 0;
    for(i = 1; i < X; i++) {
        if (getT(i,time1) > maxT && getT(i, time1) < tol) {
            maxT = getT(i, time1);
            position = i * deltaX;
        }
    }
  // printf("%lf\t%15.15lf\n",position,maxT);
    if (time1 < 200){
    //    printf("%lf\t%15.15lf\n",position*10,maxT);
    }
    waveFront[time1] = position*10;    
}

void calculateFlux(int time1){
    int i = 0;
    for (i = 0; i < X; i++) {
        F[i][time1] = - c*EF[i] * (E[i + 1][time1] - E[i][time1])/ deltaX;
           
    }
}

double convertLineToDouble(char* str, int len) {
    char delim[] = ":";
    char *ptr = strtok(str, delim);
    ptr = strtok(NULL, delim);
    double dble = atof(ptr);
    return dble;
}

int convertLineToInt(char* str, int len) {
    char delim[] = ":";
    char *ptr = strtok(str, delim);
    ptr = strtok(NULL, delim);
    int dble = atoi(ptr);
    return dble;
}

void update_dt() {
    int i, j;
    double T1, T2, max_T = 0, tmp= 0, min_T;
    double dt_tag = 0.0;
    double delta_temp = 0.0;
    //for min_T we need the max_T
    currentTime[currentTimeStep] = currentTime[currentTimeStep - 1] + deltaT;
    return ;
    for (i = 0; i < X; i++) {
        T1 = getT(i, currentTimeStep - 1);
        if ( T1 > max_T) {
            max_T = T1;
        }
    }
    min_T = max_T * 10E-2;

    for(i = 0; i < X; i++) {
        T1 = getT(i, currentTimeStep - 1);
        T2 = getT(i, currentTimeStep);
        delta_temp = fabs(T2 - T1) / (T2 + min_T);
        if (delta_temp > tmp ) {
            tmp = delta_temp;
        }
    }

    dt_tag = d_frac * (deltaT) / delta_temp;
    deltaT = Min(dt_tag, 1.1*deltaT);
    //printf("%10e\t%10e\n", delta_temp, deltaT);
}

double Min(double xx, double yy) {
    if ( xx > yy) {
        return yy;
    } else {
        return xx;
    }
}