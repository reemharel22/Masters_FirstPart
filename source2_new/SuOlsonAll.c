
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
#define NN (((X*2) + 1))
//#define NN X
//#define NN 10
float epsilon = 1e-20;
void buildABLambdaT(int XX,int NX ,float (*EF)[X], float E[X + 1][N],float F[X][N],float [X][N],int j);
void copyFromSolutionP1(float*solve,float(*mat)[X + 1][N],float (*m)[X][N],int j);
void copyFromSolutionDiff(float*solve,float(*mat)[X + 1][N],float (*m)[X][N],int j);
void sendToFileE(int);
void sendToFileT(int);
void PredictorCorrectorSolution(int times,int i,void(*f)(),void(*a)(),void(*)(),void(*)());
void constructLUDP1(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]);
void constructLUDDiff(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]);
void constructLUDP1AB(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]);
void constructLUDP1MUAB(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]);
void constructLUDDiffMUB(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]);
void CalculateT(int i,float deltaT);
void ApplyTandSourceP1(int i,float deltaX,float deltaT);
void ApplyTandSourceDiff(int i,float deltaX,float deltaT);
inline float getInitValue();
int checkConverged(int i);
float calculateMu(int space,int time1);
float calculateKappa(float);
float calculateA(float weff);
int setUpProgram(int,char*arg[]);
float getSource(int space,int time1);
float getCv(int,int);
float getT(int space,int time1);
void setUpInitialCondition();
float getTH();
float Avg1(float,float);
float Avg2(float,float);
void applyBC();
int currentTimeStep = 0;
float F[X][N],T[X + 1][N];
//@@@CHANGED TO 2*NN +1
float L[NN],U[NN],mainD[NN],E[X+1][N],EF[X+1],D[X],solve[NN],Weff[X+1][N];
float a[NN + 1], b[NN + 1], si[NN + 1], r[NN + 1];
float x0 = 0.5;
float t0 = 10.0;
float A[X+1];
float eps = 1.0;
float opacity;
float B[X+1];
float c = 3E10;
// 50 mg/cm^3
float rhoSilicon = 50;
//convert it to g/cm^3
float rho = 0.5;
float
//float c = 1;
float arad = 7.56E-15;
float P1 = 0;
float Cv = 0;
//-silicon is 3.53
float alpha = 3.53;
float beta = 1.1;
float musio = 0.09;
float s_lambda = 0.75;
float s_f = 8.78;
//in g/cm^3
float s_g = 1.0/9175.0;
//float alpha = 1.0;
int Classic = 0;
//float A = 3,B = 3;
float lambdaT;
float previousWeff = 0;
//in cm
float deltaX = 0.001;
//in nano
float deltaT = 0.01;
float TH = 1.0;
int constOpacity = 0;

int main(int argc,char *argv[]) {
  float a,b,d;
  FILE*fp;
  int k,p,h,i=0,j=0;
  void (*funcptr) (int ,int ,float(*)[X],float[X][N],float[X+1][N],float[X][N],int);
  void (*BuildLUD)(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]);
  void (*applyTandS)(int,float,float);
  void (*copySolve)(float*,float(*mat)[X+1][N],float (*m)[X][N],int j);
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

  for ( i = 0; i < N; i++) {
        for ( j = 0; j < X; j++) {
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000) {
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
              //if (2*j == 1 || 2*j == 10 || 2*j == 17 || 2*j == 31 || 2*j == 45 || 2*j == 50 || 2*j == 56 || 2*j == 75 || 2*j == 100 || 2*j == 133 || 2*j == 177)
                if (constOpacity) {
                  printf("%f\t",E[j][i]);
                }else
                 printf("%lf\t",pow(T[j][i]/arad,0.25));
                 //
                   // printf("%f\t",pow(E[j][i],0.5));
            }
        }
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000)
        printf("\n");
    }
    if (constOpacity == 1) {
        sendToFileE(p);
    }
   else {
       sendToFileT(p);
   }
  return 0;
}

//copies from solve to the matrix
void copyFromSolutionP1(float*solve,float(*E)[X + 1][N],float(*F)[X][N],int j) {
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

void copyFromSolutionDiff(float*solve,float(*E)[X + 1][N],float(*F)[X][N],int j) {
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
      fprintf(fp, "%f ",deltaT*c*(i) );
    }
    fprintf(fp, "\n");
    for ( j = 0; j < N; j++) {
      for ( i = 0; i < X + 1; i++) {
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
 } else if ( p == 11) {
    fp = fopen("../data/Temp/SuOlsonDiffusionAsymptoticData.txt","w");
  } else if ( p == 12) {
   fp = fopen("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt","w");
  }
    for ( i = 0; i < N; i++) {
            fprintf(fp, "%f ",deltaX*(i) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%f ",deltaT*c*(i) );
    }
    fprintf(fp, "\n");
    for ( j = 0; j < N; j++) {
      for ( i = 0; i < N; i++) {
            fprintf(fp,"%f ",pow(T[i][j]/arad,0.25));
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
}

void buildABLambdaT(int TTT,int nothing ,float (*EF)[X], float E1[X + 1][N],float F1[X][N],float T1[X][N],int j) {
  int i,k;
  for (i = 0; i < X + 1; i++) {
      float weff = calculateWeff(i,j);
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
    float copySolution[NN];
    for(k = 0; k < NN; k++){
        copySolution[k] = solve[k];
    }
    //we first do the basis, where we calculate E*,F*
    (*f)(X,N,&EF,E,F,T,i-1);//build EF or FL
    (*BuildLUD)(&L,&U,&mainD); // build LUD
    
    for (k = 0; k < NN; k++) {
        if (k % 2 == 1){
        //    float mu = calculateMu(k /2, currentTimeStep-1);
          //  solve[k] = solve[k] * A[k/2] ;
        } 
    }
    
    for (k = 0; k < NN; k++) {
        a[k + 1] = L[k];
        //printF()
        b[k + 1] = mainD[k];
        si[k + 1] = U[k];
        r[k + 1] = solve[k];
    }
    float *ss = malloc((NN + 2) * sizeof(float));

    tridag(a, b, si, r, ss, NN);
    for (k = 0; k < NN; k++) {
        solve[k] = ss[k + 1];
        if (solve[k] < 0 ){ 
            //solve[k] = getInitValue();
        }
    }
    free(ss);
    //now that we solved u(x,t+1), we will copy it to E.
    (*copySolve)(solve,&E,&F,i); //copy solution
    //note, when we solve the real E(n+1) we need copySolution
    CalculateT(i,deltaT);//we calculate Tn+
   //we apply to solve Tn+1 and the src for the next step  
         (*ApplyTS)(i,deltaX,deltaT);
    return;
    /*for ( j = 0; j < NN; j++) {
        solve[j] = copySolution[j];//copySolution contains E(n),F(n)
    }
     (*ApplyTS)(i,deltaX,deltaT);
    
    (*f)(X,N,&EF,E,F,T,i);//build EF or FL
    (*BuildLUD)(&L,&U,&mainD);
    for (k = 0; k < NN; k++) {
        if (k % 2 == 1){
            float mu = calculateMu(k /2, currentTimeStep-1);
            solve[k] = solve[k] * A[k/2] ;
        } 
    }
    
    for (k = 0; k < NN; k++) {
        a[k + 1] = L[k];
        b[k + 1] = mainD[k];
        si[k + 1] = U[k];
        r[k + 1] = solve[k];
    }
    ss = malloc((NN + 2) * sizeof(float));

    tridag(a, b, si, r, ss, NN);
    for (k = 0; k < NN; k++) {
        solve[k] = ss[k + 1];

    }
    free(ss);

    (*copySolve)(solve,&E,&F,i); //copy solution
    CalculateT(i,deltaT);
    (*ApplyTS)(i,deltaX,deltaT);*/

}

void constructLUDP1(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]) {
    int i,j = 0;
    float opacity;
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
          float opacityprev = getOpacity(k-1,currentTimeStep-1);
          float opacitycurr = getOpacity(k,currentTimeStep-1);
          opacity = getOpacity(k - 1,currentTimeStep-1);
          opacity = Avg2(opacitycurr,opacityprev);
          (*mainD)[i] = 3.0 + 3.0*opacity*deltaT*c;
          j++;
      }
    }
}

void constructLUDP1AB(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]) {
    int i,j = 0;
    float opacity;
    float bAvg,muAvg,aAvg;
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
          float opacityprev = getOpacity(k-1,currentTimeStep-1);
          float opacitycurr = getOpacity(k,currentTimeStep-1);
          opacity = Avg2(opacitycurr,opacityprev);
          (*mainD)[i] = A[j] + opacity*B[j]*deltaT*c;
          j++;
      }
    }

}

void constructLUDP1MUAB(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]) {
    int i,j = 0;
    constructLUDP1AB(L,U,mainD);
     //(*L)[0] = 0.0;
    j = 0;
    int k = 0;
    for ( i = 0; i < NN; i++) {
        if (i % 2 == 0) {//E
            k++;
        } else {
            float muprev,mucurrent;
            float bAvg,muAvg,aAvg;
            muprev = calculateMu2(calculateWeffNonAvg(j, 
            currentTimeStep - 1));
            mucurrent = calculateMu2(calculateWeffNonAvg(j , 
            currentTimeStep - 1));
            float udiv = mucurrent;
            float ldiv1 =  muprev;
            
            (*U)[i] = (c*udiv*c*deltaT)/(deltaX);
            (*L)[i] = -(c*ldiv1*c*deltaT)/(deltaX);
            float opacityprev = getOpacity(k-1,currentTimeStep-1);
            float opacitycurr = getOpacity(k,currentTimeStep-1);
            opacity = Avg2(opacitycurr,opacityprev);
            float mu = calculateMu(j , currentTimeStep - 1);
            (*mainD)[i] = A[j]*mu + mu*opacity*B[j]*(c*deltaT);
            if (currentTimeStep == 273) {
               // printf("%lf\t%lf\t%lf\t%d\n",udiv,ldiv1,mu,i/2);
               
            }
            j++;
        }
    }
    (*L)[0] = 0.0;

}

void constructLUDDiff(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]) {
    int i,j = 0;
    float opacity;
    float lambda = deltaT/(deltaX*deltaX);
    for (i = 0; i < NN-1; i++) {
        if (i != X-1) {
           (*U)[i] = -lambda *c*((EF[i] + EF[i+1])/2.0);
      }
    }

    (*L)[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        if (i != 0 ) {
            (*L)[i] = -lambda*c *( (EF[i-1]+EF[i])/2.0);
        }
    }
    j = 0;
    for (i = 0; i < NN; i++) {
        opacity = getOpacity(i,currentTimeStep-1);
        (*mainD)[i] = 1 + deltaT*c*opacity + lambda*c*((EF[i+1]+2*EF[i] + EF[i-1])/(2.0));
    }
    //@@@ added 2*deltaX
    float bb = 0;
    if (constOpacity == 0 ) {
        bb = deltaX;
    }
    (*mainD)[0] = 1 + deltaT*c*getOpacity(0,currentTimeStep-1) + lambda*c*(bb + EF[0] + EF[1])/2.0;
  //  printf("%lf\t",EF[1]);
    i = NN - 1;
    opacity = getOpacity(i,currentTimeStep-1);
    (*mainD)[NN-1] = 1 + deltaT*c*opacity + lambda*c*((EF[i]+2*EF[i] + EF[i-1])/(2.0));
}

void constructLUDDiffMUB(float (*L)[NN],float (*U)[NN],float (*mainD)[NN]) {
    int i,j = 0;
    float lambda = deltaT/(deltaX*deltaX);
    float mu,opacity;

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
    float bb= 0;
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

void CalculateT(int i,float deltaT) {
    int j;
    for ( j = 0; j < X + 1; j++) {
        float t = getT(j,i-1);
        float cap = getCv(j,i-1);
        float coeff = (getOpacity(j,i-1) * 4.0  * pow(t,3) * arad)
        / (cap);
        T[j][i] = ((T[j][i-1]) + deltaT*c*coeff*(E[j][i])) / (coeff*c*deltaT + 1.0);
        if ( currentTimeStep == 365) {
          //  printf("%15.15lf\t%15.15lf\t%d\t%15.15lf\n",T[j][i],E[j][i],j,calculateMu(j,i));
           // printf("%15.15lf\t%15.15lf\t%15.15lf\n\n",L[j],U[j],mainD[j]);
        }
        if (T[j][i] < 0 ){ 
           // printf("%15.15lf\t%15.15lf\t%d\n",T[j][i-1],E[j][i],j);

        }
    }  
}

void ApplyTandSourceP1(int i,float deltaX,float deltaT) {
    int j;
    int k = 0;
    float Src;
    for ( j = 0; j < X  ; j++) {
        Src = getSource(j,i);
        solve[2*j] += getOpacity(j,i)*deltaT*c*T[j][i]+ Src*deltaT*c;
    }
    //@@@added this
    if (!constOpacity) {
      applyBC();
    }
}

void ApplyTandSourceDiff(int i,float deltaX,float deltaT) {
    int j;
    int k = 0;
    float Src;
    for ( j = 0; j < X; j++) {
        Src = getSource(j,i);
        solve[j] += getOpacity(j,i)*deltaT*c*T[j][i]+ Src*c*deltaT;
    }
    if (!constOpacity) {
        float l = deltaT / (deltaX*deltaX);
        solve[0] += (l*2.0*c*deltaX*getFinc())/(c);
    }
}

int checkConverged(int j) {
    int i;
    float resuSumE = 0.0;
    float resuSumF = 0.0;
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

float getdx() {
    return deltaX;
}

float getdt() {
    return deltaT;
}

float calculateWeff(int space,int time1) {
    float weff;
    float tt,ee;
    float opacity = getOpacity(space,time1);
    float Src = getSource(space,time1);
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
        float op1 = getOpacity(space-1 ,time1);
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

float calculateWeffNonAvg(int space,int time1) {
    float weff;
    float tt,ee;
    float opacity = getOpacity(space,time1);
    float Src = getSource(space,time1);
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

float calculateKappa(float weff) {
    if (weff < 0.01) {
        return 1.0;
    } else if (0.01 <= weff && weff <= 0.45) {
        float aa = exp(-2.0/(weff));
        float bb = (24.0 + 20.0*weff + 3.0*pow(weff,2.0)) / (pow(weff,2.0) );
        float xyz = (1.0- (4.0*aa*(1.0 + aa*((4.0 - 2.0*weff)/(weff)) + bb*exp(-4.0/(weff)) ) ));
        return xyz;
    } else if (0.45 < weff && weff < 1.0) {
        return (1.0-weff) * calculateB(weff);
    } else  {
        return (weff-1.0) * calculateB(weff);
    }
}

float calculateB(float weff) {
    //float b1;
    
    if ( 0.59 <= weff && weff <=0.61) {
     return 1.0 / (0.80054 - 0.523*weff);
    } else {
        float xyz = (1.0 + weff) / 0.40528473;
        float aa = (0.1326495 + weff*(0.03424169 + weff*(0.1774006 - weff)));
        float ba = (0.3267567 + weff*(0.1587312 - weff*(0.5665676 + weff)));
        return (aa*xyz)/ba;
    }
}

float calculateA(float weff) {
    if (0.55 <= weff && weff <= 0.65) {
      return 0.96835 - 0.437*weff;
    } else {
        float aa = ( 0.247 * (0.433 + 0.421*weff -2.681*weff*weff
             - 1.82*pow(weff,3.0) + 4.9*pow(weff,4.0) -1.058*pow(weff,5.0)
              + 2.56*pow(weff,6.0) ) );
        float ba =  pow(0.327 + 0.159*weff - 0.567 * pow(weff,2.0) - pow(weff,3.0)  ,2.0);
      return aa/ba ;
    }
}

float calculateMu(int space,int time1) {
    float weff = calculateWeff(space,time1);
    return calculateMu2(weff);
    float kappa = calculateKappa(weff);
    if (weff < 0.01) {
        return 1.0;
    } else if(0.01 <= weff && weff < 0.999) {
        if (kappa > 0) {
            float a = -weff/(2.0*kappa);
            float b = a*(log(1.0-kappa + 1e-15));
            return b;
        } else {
            return 1;
        }
    } else if (0.999<= weff && weff <= 1.001) {
        float a = 8.3548 + 1.5708 + weff;
        float b = 2.1228 + 2.4674*weff;
        return log(a/b);
    } else {
        float a = weff/(2.0*kappa + 1e-15);
        return a*log(1.0+kappa);
    }
}

float calculateMu2(float weff) {
    float kappa = calculateKappa(weff);
    if (weff < 0.01) {
        return 1;
    } else if(0.01 <= weff && weff < 0.999) {
        if (kappa > 0) {
            float a = -weff/(2.0*kappa);
            float b = a*(log(1.0 - kappa + 1e-15));
            return b;
        } else {
            return 1.0;
        }
    } else if (0.999<= weff && weff <= 1.001) {
        float a = 8.3548 + 1.5708 + weff;
        float b = 2.1228 + 2.4674*weff;
        return log(a/b);
    } else {
        float a = weff/(2.0*kappa /*+ 1e-15*/);
        return a*log(1.0 + kappa);
    }
}

float getOpacity(int space,int time1) {
    if (constOpacity) {
        return 1.0;
    } else if (constOpacity == 0){
        if (space >= X) {
          space = X-1;
        }
        float v = T[space][time1];
        float t  = pow(v/arad,0.25);
        
        if (v < epsilon || t < epsilon) {
            t = 0.1;
            //printf("hm");
        }
        float a = 1.0/(pow(t/(getTH()),3));
        if(currentTimeStep == 854 || a > 1000) {
          //  printf("%25.25lf\t%lf\n",t,a);
        }
        return a;
    }else {
        double xx = s_g * pow(getT(space,time), alpha) * pow(rho, -s_lambda);
        return rho / xx; 
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

float getSource(int space,int time1) {
    //we are in a const opacity i.e src is 1
    if (constOpacity) {
         if (! (space*deltaX >= x0 || time1*c*deltaT >= t0))
        {
            return 1.0;
        }
        return 0.0;
    } else {
        return 0;
    }
}

void setUpInitialCondition() {
    float Src;
    int i,j;
    deltaT = deltaT/c;
    TH = 1.0;
    t0 = t0;
    //alpha = 4.0*arad/eps;
    for ( i = 0; i < NN; i++) {
        solve[i] = 0;
    }

    if (P1) {
      for ( i = 0; i < X + 1; i++) {
        Src = getSource(i,currentTimeStep);
        solve[2*i + 1] = Src*deltaT*c;
      }
    }
    else {//diffusion
      for ( i = 0; i < X; i++) {
           Src = getSource(i,currentTimeStep);
          //@@@ add if & else
          if (i == 0 && constOpacity == 0) {
            float l = deltaT/(deltaX*deltaX);
                    solve[0] += (l*2.0*c*deltaX*getFinc())/(c);
          } else {
              solve[i] = Src*deltaT*c;
          }
    }
  }
  //#pragma omp parallel for default(shared)
    if (constOpacity) {
        for (i = 0; i < X; i++) {
          for ( j = 0; j < N; j++) {
              T[i][j] = E[i][j] = 0;
          }
        }
    } else {
        for (i = 0; i < X + 1; i++) {
          for ( j = 0; j < N; j++) {
           // float tm = pow(10,-1.25)*5*getTH();
             // T[i][j] = E[i][j] = arad*pow(tm,4);
                T[i][j] = 190ev;
                E[i][j] = getInitValue();

          }
        }
    }
    for (i = 0; i < NN; i ++) {
        solve[i] = E[0][0];
    }
    for (i = 0; i < X; i++) {
      for ( j = 0; j < N; j++) {
          F[j][i] =  getInitValue();
      }
    }
    for ( i = 0; i < X; i++) {
        D[i] = EF[i] = (float)1.0/(3.0);
    }
    for ( i = 0; i < X; i++) {
      A[i] = B[i] = 3.0;
    }
    if (!constOpacity) {
        currentTimeStep = 1;
        applyBC();
      }
}

float getCv(int space, int time1) {
    if (constOpacity) {
        float a = getT(space,time1);
        return (alpha*pow(a,3));
    } else {
        return 4.0*arad*pow(getTH(),3);
    }
}

float getT(int space,int time1) {
    return pow(T[space][time1]/arad,0.25);
}

float getTH(){
    return 1.0;
}

float getFinc(){
    return (c*arad*pow(getTH() * 5,4)/4.0);
}

void applyBC() {
  if (!constOpacity) {
    float mu = calculateMu2(calculateWeff(0,currentTimeStep ));
    solve[0] += (2.0*getFinc() - c*solve[0]*mu)*(deltaT/deltaX);
  } 
}

float Avg1(float xx,float yy) {
    return (xx + yy) / 2.0;
}

float Avg2(float xx,float yy) {
    if ( xx == 0.0 || yy == 0.0) {
        return Avg1(xx,yy);
    }
   return (2.0 * xx * yy)/(xx + yy);
}
inline float getInitValue() {
    if (constOpacity == -1) {
        //todo ev...
        return arad * pow(190,4)
    } else {
      float tm = pow(10,-1.25)*5*getTH();
      return arad*pow(tm,4);
      }
}