#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "SuOlsonAll.h"
#include <malloc.h>
#include <string.h>
#include "utils/nrutil.h"
#include "utils/nrutil.c"
#include "utils/tridag.c"
#include "tridFunc.h"
#define N 1500 // isn't relevant anymore
#define X 1500 // space
#define NN 1500

double epsilon = 1e-20;
double get_time();
void sendToFileE(int);
void sendToFileT(int);
double Max(double, double);
void PredictorCorrectorSolution(int times,int i,void(*f)(),void(*a)(),void(*)());
void constructLUDDiff();
void sendToFileW(int p);
void CalculateT();
void ApplyTandSourceDiff(int i,double deltaX,double deltaT);
inline double getInitValue();
int setUpProgram(int,char*arg[]);
double getSource(int space,int time1);
double getT(int space,int time1);
void setUpInitialCondition();
double getTH();
double Avg1(double,double);
void diagnostics();
double convertLineToDouble(char* str, int len);
int convertLineToInt(char* str, int len);
double* convertLineToArray(char* str, int len);
void update_opacity();
double Avg2(double,double);
void update_dt();
void applyBC(int);
double Min(double, double);
int currentTimeStep = 0;
double F[X][N],T[X + 1][N];
double prev_time;
//@@@CHANGED TO 2*NN +1
double L[NN + 1],U[NN + 1],mainD[NN + 1],E[X],EF[X],D[X],solve[2 * NN + 1],Weff[X+1];
double opac[X];
double E_old[X], E_current[X], V_old[X], V_current[X]; // V ofc is the mat tempretaure and e is the rad temp
double waveFront[X];
double currentTime;
int problem;
double *diag;
int num_diagnostics;
double initV = 1E-50;
double TH[2][759];
double a_[NN + 1], b_[NN + 1], si_[NN + 1], r_[NN + 1];
double x0 = 0.5;
double t0 = 10.0;
double sig_factor = 1;
double eps = 1.0;
double opacity;
double sigma_boltzman = 5.670373e-5;
double c = 3E10;
int BC = 0;
double rho = 0.05;
double arad = 7.56E-15;
double d_frac = 0.1;
double alpha;
double beta ;
double mu_sio = 0.1;
double s_lambda = 0.75;
double s_f = 8.8E13;
double s_g = 1.0/9175.0;
int Classic = 0;
double lambdaT;
double bb = 0;
double t_s = 0, t_tot = 0;
double previousWeff = 0;
double time_finish;
double deltaT = 0.01;
double deltaX = 0.01;
double th;
int constOpacity = 0;

int main(int argc,char *argv[]) {
    double a,b,d;
    FILE*fp;
    int k,p,h,i=0,j=0;
    omp_set_num_threads(4);
    void (*funcptr) (int  ,double(*)[X], double [X]);
    void (*BuildLUD)();
    void (*applyTandS)(int,double,double);
    double t1 = get_time();
    // now we have a trid-matrix size 2N on 2N
    //ROWS IS FOR Space ie i const, j not. you move in T axis
    //COLS IS FOR TEperture ie j const, i not. you move in X axis
    p = setUpProgram(argc,argv);//read datafile
    
    funcptr = &buildDDiff; // funcptr is how we build the diffusion coefficient.
    BuildLUD = &constructLUDDiff;// build the tridiagonal matrix
    applyTandS = &ApplyTandSourceDiff;//Apply T and source to the equation of energy..

    setUpInitialCondition();
    currentTimeStep = 1;

    double kk = c;
    if (problem == 2) {
        kk = 1E9;
    }
    while (currentTime * kk < time_finish) {
        if (currentTimeStep % 1000 == 0) {
            printf("Current time: %10e\t Time step: %d\n", currentTime, currentTimeStep);
        }
        PredictorCorrectorSolution(1, currentTimeStep, funcptr,BuildLUD,applyTandS);

        currentTimeStep++;
    }
    printf("Time took: %10e\n",get_time() - t1);
    printf("Time took for wheel change: %10e\n",t_tot);
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

/**
 * This is the timestep of the program and it goes as the following:
 * 1. update_opacity: Calculate the opacities of this time step via E_Old and T_Old aka E(i,n-1)
 * 2. (*f): Calculate the diffusion coefficient - happends in tridFunc.c .
 * 3. ApplyTandSourceDiff: As the energy equation is with T and sources, we will want to add them to
 * the solution vector, this is done here. (Ax=b, we add these terms to the b which is called solve[] here).
 * 4. applyBC: We apply the boundary condition to the equation, in many cases we apply it to the solve array which is the b vector
 * 5. (*BuildLUD): We build the tridiagonal matrix.
 * 6. solveTriagonal: We solve the Ax=b equation, the solution is stored in the b vector (solve array).
 * 7. We copy the solution to E_Current
 * 8. We calculate the temperature.
 * 9. diagnostics: To see if we need to write to a file..
 * 10. update_time: calculate new dt
 * 
 * **/
void PredictorCorrectorSolution(int times, int i, void(*f)(),void(*BuildLUD)(),void(*ApplyTS)()) {
    int j,k=0,p=0;
    t_s = get_time();

    update_opacity();
    (*f)(X,&EF, opac); //build EF or FL
    
    t_tot += get_time() - t_s;

    ApplyTandSourceDiff(i - 1, deltaX, deltaT);

    applyBC(i - 1);
    (*BuildLUD)(&L,&U,&mainD); // build LUD

    solveTriagonal(NN, &solve, L, U, mainD);

    for (j = 0; j < X; j++) {
        E_current[j] = solve[j];
    }
    CalculateT();

    diagnostics();
    update_dt();
    for ( j = 0; j < X; j ++) {
        E_old[j] = E_current[j];
        V_old[j] = V_current[j];
    }

    return;

}

//Calculate the new opacities via getOpacity.
void update_opacity() {
    int i = 0;
    for (i = 0 ; i < X; i++) {
        opac[i] = getOpacity(i, 0);
    }
}

//Calculates the opacity with the equation sigma = kappa * rho.
double getOpacity(int space, int time1) {
    double t_galpha = (pow(rho, s_lambda + 1) 
                        / (s_g * pow(getT(space, time1), alpha)));
    return t_galpha;
}

//Applies the temperature and source to the energy equation
void ApplyTandSourceDiff(int i,double deltaX,double deltaT) {
    int j;
    int k = 0;
    double Src = 0;
    for ( i = 0; i < X; i++) {
        Src = getSource(j, 0);
       // solve[i] += sig_factor * getOpacity(i, 0)* deltaT *c*V_old[i] + Src*c*deltaT;
        solve[i] += sig_factor * opac[i] * deltaT *c*V_old[i] + Src*c*deltaT;
    }
    //applyBC(currentTimeStep - 1);
}

double getSource(int space,int time1) {
    //we are in a const opacity i.e src is 1
         if ( (space*deltaX < x0 && currentTime * c < t0)) {
            return 1.0;
        }
        return 0.0;
}

//applies the bc to the equation
void applyBC(int time1) {
    if (constOpacity == 0) {
        //double mu = calculateMu2(calculateWeff(0,currentTimeStep ));
        //solve[0] += getFinc() * deltaT * 2.0/(c * deltaX);
    double l = deltaT / (deltaX*deltaX);
        solve[0] += l*c*arad*deltaX/2.0;
        //solve[0] += (2.0*getFinc() - c*solve[0]*mu)*(deltaT/deltaX);
    } else if (constOpacity == -1) {
        if (BC == 1) {
            solve[0] += ( (2.0*getFinc() * EF[0] * (deltaT/deltaX)) / (EF[0] + deltaX*0.5)); // bc avner 2
        } else if ( BC == 2) {
            //solve[0] += arad * c * pow(getTH(time1) ,4) * deltaT/(2.0*deltaX);
            solve[0] += 2.0*getFinc() * deltaT / deltaX;  // bc avner 2
        } else if ( BC == 3) {
            solve[0] = 2.0*getFinc() * deltaX/(c*EF[0]);
        } else if ( BC == 4) {
            solve[0] = arad*pow(getTH(time1), 4)/(1.0); 
        }else if ( BC == 5) {
        //    double coeff = 100*11605*pow(currentTime[time1]/1E-9, 0.1408);
            double coeff = 100 * 11605 * pow(currentTime / 1E-9, 0.1163);
            coeff = Max(300, coeff);
            solve[0] = arad * pow(coeff,4);
            
        }
    } else if (constOpacity == 1) {
        // solve[0] += getFinc() * deltaT * 2.0/(c * deltaX);
    }
}

//Builds the tridiagonal A matrix
void constructLUDDiff() {
   int i,j = 0;
    double opacity;

    double lambda = deltaT/(deltaX*deltaX);
    for (i = 0; i < NN-1; i++) {
        if (i != X-1) {
            U[i] = -lambda *c* (EF[i + 1] + EF[i])/2.0;
        }
    }

    L[0] = 0.0;
    j = 1;
    for ( i = 1; i < NN; i++) {
        L[i] = -lambda * c * (EF[i -1] + EF[i])/2.0;
    }
    j = 0;
    for (i = 0; i < NN; i++) {
        //opacity = sig_factor * getOpacity(i, 0);
        opacity = sig_factor * opac[i];
        mainD[i] = 1 + deltaT*c*opacity + lambda*c*(EF[1 + i] + 2*EF[i] + EF[i - 1])/(2.0);
    }
    //@@@ added 2*deltaX
    //    double bb = 0;
    // doesn't work for su olson.
    if (BC == 1) {
        bb = -2*pow(EF[0], 2) / (EF[0] + 0.5*deltaX);// avner bc 1
        bb += 2.0*EF[0]; 
    }
    mainD[0] = 1.0 + deltaT*c*sig_factor * getOpacity(0, 0) + lambda*c*(bb + EF[0] + EF[1])/2.0;
    i = NN - 1;
    if (BC == 3) {
        mainD[0] = 1;
        U[0] = deltaX/(2.0*EF[0]) - 1 ;
    } else if ( BC == 4 || BC == 5) {
        mainD[0] = 1;
        U[0] = 0;
    }
   // printf("%lf\t",EF[0] );
    opacity = sig_factor * getOpacity(i,0);
    mainD[NN-1] = 1.0 + deltaT*c*opacity + lambda* c *(EF[i]+2*EF[i] + EF[i-1])/(2.0);
}

//Calculates the temeprature
void CalculateT() {
    int i;
    int stop = 3;
    double coeff;
    for ( i = 0; i < X; i++) {
        double t = getT(i, 0);
        double cap = getCv(i, 0);
       // double coeff = (sig_factor * getOpacity(i, 0) * 4.0  * pow(t, 3) * arad)
       // / (cap);
       coeff = (sig_factor * opac[i] * 4.0  * pow(t, 3) * arad)
        / (cap);
        //coeff = deltaT*c*coeff;

        V_current[i] = ((V_old[i]) + (coeff*c*deltaT*E_current[i])) /
                                    (1.0 + deltaT*c*coeff);
    }  
}

//calculates heat capacity via rho*(de/dT)
double getCv(int space, int time1) {
    double abbcdef = beta * s_f * pow(getT(space, time1), beta - 1 ) * pow(rho, -mu_sio + 1.0);
    return abbcdef;
}

double getdx() {
    return deltaX;

}

double getdt() {
    return deltaT;
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
        } else if (strstr(line, "x0:") != NULL) {
            x0 = convertLineToDouble(line, len);
        } else if(strstr(line, "t0:") != NULL) {
            t0 = convertLineToDouble(line, len);
        } else if(strstr(line, "g:") != NULL) {
            s_g = convertLineToDouble(line, len);
        } else if(strstr(line, "f:") != NULL) {
            s_f = convertLineToDouble(line, len);
        } else if(strstr(line, "mu:") != NULL) {
            mu_sio = convertLineToDouble(line, len);
        } else if(strstr(line, "lambda:") != NULL) {
            s_lambda = convertLineToDouble(line, len);
        } else if(strstr(line, "deltaT:") != NULL) {
            deltaT = convertLineToDouble(line, len);
        } else if(strstr(line, "deltaX:") != NULL) {
            deltaX = convertLineToDouble(line, len);
        } else if(strstr(line, "alpha:") != NULL) {
            alpha = convertLineToDouble(line, len);
        } else if(strstr(line, "beta:") != NULL) {
            beta = convertLineToDouble(line, len);
        } else if(strstr(line, "rho:") != NULL) {
            rho = convertLineToDouble(line, len);
        } else if(strstr(line, "TH:") != NULL) {
            th = convertLineToDouble(line, len);
        } else if(strstr(line, "Temp_0:") != NULL) {
            initV = convertLineToDouble(line, len);
        } else if(strstr(line, "dfrac:") != NULL) {
            d_frac = convertLineToDouble(line, len);
        } else if(strstr(line, "BC:") != NULL) {
            BC = convertLineToInt(line, len);
        } else if(strstr(line, "sig_factor:") != NULL) {
            sig_factor = convertLineToDouble(line, len);
        } else if(strstr(line, "diagnostics:") != NULL) {
            diag = convertLineToArray(line, len);
        } else if(strstr(line, "timefinish:") != NULL) {
            time_finish = convertLineToDouble(line, len);
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
                continue;
            } else if (strstr(line,"AB") != NULL) {
                p = 4;
                continue;
            } else if (strstr(line,"EF") != NULL) {
                if (strstr(line,"Kershaw") != NULL) {
                    p = 0;
                    continue;
                } else if (strstr(line,"LP") != NULL) {
                    p = 1;
                    continue;
                } else if (strstr(line,"Minerbo") != NULL) {
                    p = 2;
                    continue;
                }
            }
            p = 3;
        }
    }
    if (problem == 0) {
        deltaT = deltaT / c; //units
        bb = 0;
    } else if (problem == 2) {
      //  deltaT = deltaT/c;
        deltaT = deltaT * 1E-9; //units
        bb = deltaX;
        s_f = s_f / (pow(1160500.0, beta));
        
        s_g = s_g /  pow(1160500.0, alpha);
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

void setUpInitialCondition() {
    double Src;
    int i,j;

    //initial value for the energy and temprature
    for ( j = 0; j < N; j++) {
        E_old[j] = V_old[j] = pow(initV, 4) * arad;
    }

    for ( i = 0; i < X; i++) {
        D[i] = EF[i] = (double)1.0/(3.0 * getOpacity(i, 0));
        solve[i] = E_old[i];
    }

    char buff[255];

    currentTime = deltaT;
    prev_time = 0;
    FILE*fp;
    if (problem == 0) {
        fp = fopen("../data/SuOlsonTemp.txt","w");
    } else if (problem == 1) {
        fp = fopen("../data/OlsonTemp.txt","w");
    } else if (problem == 2) {
        fp = fopen("../data/BackTemp.txt" ,"w");
    }
    for (i = 0; i < X; i ++) {
        fprintf(fp, "%10e ",deltaX * i);
    }
    fprintf(fp, "\n");
    fclose(fp);
     if (problem == 0) {
        fp = fopen("../data/SuOlsonEnergy.txt","w");
    } else if (problem == 1) {
        fp = fopen("../data/OlsonEnergy.txt","w");
    } else if (problem == 2) {
        fp = fopen("../data/BackEnergy.txt","w");
    }
     for (i = 0; i < X; i ++) {
        fprintf(fp, "%10e ",deltaX * i);

    }
    fprintf(fp, "\n");
    fclose(fp);
}

// GETS THE REAL TEMPERATURE, Time1 = 0 we return previous temprature, time1 = 1 we reutrn the current one
double getT(int space,int time1) {
    if (time1 == 0) { // we will return the old temp
        return pow(V_old[space] / arad,0.25);
    } else {
        return pow(V_current[space] / arad,0.25);
    }
    
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


void update_dt() {
    int i, j;
    double T1, T2, max_T = 0, tmp= 0, min_T;
    double dt_tag = 0.0;
    double delta_temp = 0.0;
    //for min_T we need the max_T
    prev_time = currentTime;
    currentTime = currentTime + deltaT;
    for (i = 1; i < X; i++) {
        T1 = getT(i, 0);
        if ( T1 > max_T) {
            max_T = T1;
        }
    }
    if (problem == 2) {
        max_T = 190 * 11605;
    }
   // max_T *= 10;
    //max_T = 100 * 11605;
    min_T = max_T * 1E-3;

    for(i = 1; i < X; i++) {
        T1 = getT(i, 0 );
        T2 = getT(i, 1);
        delta_temp = fabs(T2 - T1) / (T2 + min_T);
        if (delta_temp > tmp ) {
            tmp = delta_temp;
        }
    }

    dt_tag = d_frac * (deltaT) / tmp;
    deltaT = Min(dt_tag, 1.1*deltaT);
   // if (deltaT > 1E-11) {
    //    deltaT = 1E-11;
    //}
 //   printf("dt tag: %10e\tdt: %10e\tt: %10e\tdelta Temp: %10e\n",dt_tag,deltaT, currentTime, tmp);
    if (currentTimeStep == 300) {
       // exit(1);
    }
}

double Min(double xx, double yy) {
    if ( xx > yy) {
        return yy;
    } else {
        return xx;
    }
}

double Max(double xx, double yy) {
    if ( xx < yy) {
        return yy;
    } else {
        return xx;
    }
}

void diagnostics() {
    int i;
    double kk = c;
    if (problem == 2) {
        kk = 1E9;
    }
    if (problem == 0) { // su olson we want times of 
        if (prev_time * kk <= 3.16 && currentTime * kk >= 3.16) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        } else if (prev_time * kk <= 10 && currentTime * kk >= 10) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        }
    } else if(problem == 1) { // olson 
        if (prev_time * kk <= 3 && currentTime * kk >= 3) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        } else if (prev_time * kk <= 10 && currentTime * kk >= 10) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        }
    } else if( problem == 2) {
        if (prev_time * kk <= 1 && currentTime * kk >= 1) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        } 
    }
    for ( i = 0; i < num_diagnostics; i ++) {
        if (prev_time * kk <= diag[i] && currentTime * kk >= diag[i]) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        } else if (diag[i] == 9999 && currentTimeStep % 10 == 0) {
            sendToFileE(1);
            sendToFileT(1);
          //  findWavefront(1);
        }
    }
}

void sendToFileE(int p) {
    int i = 0,j;
    FILE*fp;
    if (problem == 0) {
        fp = fopen("../data/SuOlsonEnergy.txt","a");
    } else if (problem == 1) {
        fp = fopen("../data/OlsonEnergy.txt","a");
    } else if (problem == 2) {
        fp = fopen("../data/BackEnergy.txt","a");
    }
    double kk = c;
    if (problem == 2){
        kk = 1E9;
    }
    fprintf(fp, "%d ", currentTimeStep - 1 );
    fprintf(fp, "%10e ", (prev_time) * kk );
    for ( i = 0; i < X; i++) {
            fprintf(fp,"%10e ",  pow(E_old[i] / arad, 0.25));
    }
    fprintf(fp,"\n");
    fprintf(fp, "%d ", currentTimeStep );
    fprintf(fp, "%10e ", (currentTime) * kk );
    for ( i = 0; i < X; i++) {
            fprintf(fp,"%10e ",  pow(E_current[i] / arad, 0.25));
    }
    fprintf(fp,"\n");
    printf("Diagnostic: Time: %10e\t Step: %d, %d\n", currentTime * kk, currentTimeStep);
    fclose(fp);
}

void sendToFileT(int p) {
    int i = 0,j;
    FILE*fp;
    if (problem == 0) {
        fp = fopen("../data/SuOlsonTemp.txt","a");
    } else if (problem == 1) {
        fp = fopen("../data/OlsonTemp.txt","a");
    } else if (problem == 2) {
        fp = fopen("../data/BackTemp.txt","a");
    }
    double kk = c;
    if (problem == 2){
        kk = 1E9;
    }

    fprintf(fp, "%d ", currentTimeStep - 1 );
    fprintf(fp, "%10e ",(prev_time) * kk );
    if (problem != 2){
        for ( i = 0; i < X; i++) {
            fprintf(fp,"%10e ",(getT(i, 0)));
        }
    } else {
        for ( i = 0; i < X; i++) {
                fprintf(fp,"%10e ",(getT(i, 0)/11605.0));
        }
    }
   // printf("\n");
    fprintf(fp,"\n");
    fprintf(fp, "%d ", currentTimeStep );
    fprintf(fp, "%10e ",currentTime * kk );
    if (problem != 2){
        for ( i = 0; i < X; i++) {
            fprintf(fp,"%10e ",(getT(i, 1)));
        }
    } else {
        for ( i = 0; i < X; i++) {
                fprintf(fp,"%10e ",(getT(i, 1)/11605.0));
        }
    }

    fprintf(fp,"\n");
    fclose(fp);
}

void sendToFileW(int p) {
    int i = 0,j;
    FILE*fp;
   if ( p == 8) {
      fp = fopen("../data/Temp/Back_1500_WaveFront.txt","w");
   }
    for ( i = 0; i < X; i++) {
        fprintf(fp, "%10e ",deltaX*(i) );
    }
    fprintf(fp, "\n");
    for ( i = 0; i < N; i++) {
      fprintf(fp, "%10e ",currentTime * 1E9 );
    }
    fprintf(fp, "\n");
    fclose(fp);
}

double get_time() {
  return omp_get_wtime();
}

//for datafile
double convertLineToDouble(char* str, int len) {
    char delim[] = ":";
    char *ptr = strtok(str, delim);
    ptr = strtok(NULL, delim);
    double dble = atof(ptr);
    return dble;
}

//for datafile
int convertLineToInt(char* str, int len) {
    char delim[] = ":";
    char *ptr = strtok(str, delim);
    ptr = strtok(NULL, delim);
    int dble = atoi(ptr);
    return dble;
}

//for datafile
double* convertLineToArray(char* str, int len) { 
    char delim[] = ":";
    len = -1;
    char * tmp[50];
    strcpy(tmp, str);
    char *ptr = strtok(tmp, delim);
    while (ptr != NULL) {
        len ++;
        ptr = strtok (NULL, " ");
    }
    double * arr = (double*) malloc(sizeof(double) * len);
    char *ptr1 = strtok(str, delim);
    int i = 0;
    for ( i = 0; i < len; i++) {
        ptr1 = strtok (NULL, " ");
        arr[i] = atof(ptr1);
    }
    num_diagnostics = len;
    return arr;
}