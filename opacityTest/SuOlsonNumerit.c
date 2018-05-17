#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>
#define N 10
#define X 10
#define n (N*2 + 1)
void matrixDotVector(double mat[][N],double (*vec)[N][N],int j);
void printMatrix(double mat[][N]);
void printMatrix2(double mat[][n]);
void setMatrixForE(double (*mat)[n][n],double , double,double);
void solveTriagonal(double(*solve)[n],double L[n],double U[n],double mainD[n]);
void constructLUD(double (*L)[n],double (*U)[n],double (*mainD)[n],double [n][n]);
void copyFromSolution(double*solve,double(*mat)[X+1][N],double (*m)[X][N],int j);
void copyToSolution(double(*solve)[N],double mat[N][N],int j);
double integrateX(double(*U)[N][N],int i,double);
void sendToFile(double E[N][N],double T[N][N],double,double);
double getOpacity(int space,int time1);
double getCv(int space, int time1) ;
double getFinc() ;
double getT(int space,int time1);
//void mallocMatrices(double**)
double E[X+1][N],T[X][N],matrix[n][n];
double L[n],U[n],mainD[n],solve[n],F[X][N];
int currentTimeStep = 0;
double Cv = 0;
double alpha = 1;
double arad = 7.56E-15;
double eps = 1;
double constOpacity = 1;
double c = 3E10;
int main() {
  int k,p,h,i=0,j=0;
  double sigma;
  double lambdaE,lambdaT,Tp,x0,t0;
  double b,d,Src;
  double deltaX,deltaT;
  FILE*fp;
  alpha = 4.0*arad;
  //void *()
  // now we have a trid-matrix size 2N on 2N
  //ROWS IS FOR Space ie i const, j not. you move in T axis
  //COLS IS FOR TEperture ie j const, i not. you move in X axis
  sigma = 1;
  //boundry of the source
  x0 = 0.5;
  t0 = 10.0;
  deltaX = (double)0.01;
  deltaT = (double)0.01/c;
  lambdaE = deltaX;
  lambdaT = deltaT;
  Src = 0;
  //setting up the matrices
  for ( i = 0; i < n; i++) {
      if (i%2 == 0) {
          Src = 1;
          if (i*deltaX/2 >= x0)
          {
              Src = 0;
          }
          solve[i] = Src*deltaT*c;
      } else {
          solve[i] = 0;
      }
  }
  //solve[0] = 0.5;
  for (i = 0; i < X+1; i++) {
    for ( j = 0; j < N; j++) {
        F[i][j] = T[i][j] = E[i][j] = pow(10,-5);
    }
  }
  for (i = 0; i < n; i++) {
    for ( j = 0; j < n; j++) {
        matrix[i][j] = 0.0;
    }
  }

  //sets the matrix - Au(x,t) = u(x,t-1)
  setMatrixForE(&matrix,deltaX,deltaT,lambdaT);
  //constructs the upper,lower diagonals.
  constructLUD(&L,&U,&mainD,matrix);
   
  for (i = 1; i < N; i++) {
    currentTimeStep = i;
    //if we want to change the deltaT need TODO:
    //update lambdaE & lambdaT and call seTatrixforE & constructLUD again
    solveTriagonal(&solve,L,U,mainD);
    //now that we solved u(x,t+1), we will copy it to E.
    copyFromSolution(solve,&E,&F,i);
    setMatrixForE(&matrix,deltaX,deltaT,lambdaT);
    constructLUD(&L,&U,&mainD,matrix);
    //#pragma omp parallel for default(shared)
   
    for ( j = 0; j < N; j++) {
        double ttt = getT(j,i-1);
        double cap = getCv(j,i-1);
        double coeff = (getOpacity(j,i-1) * 4.0  * pow(ttt,3) * arad) 
        / (cap);
       
      //this is where we calculate the next Temperture !
       T[j][i] = ((T[j][i-1]) + deltaT*c*coeff*E[j][i]) / (coeff*c*deltaT + 1.0);
       if (i > 0 && i < 50) {
        // printf("%lf\n", E[j][i]);
       }
     }
     //#pragma omp parallel for default(shared) private (Src)
     for ( j = 0; j < n; j+=2) {
         Src = 1;
         if (j*deltaX/2 >= x0 || i*deltaT*c >= t0)
         {
             Src = 0;
         }
          solve[j] += getOpacity(j/2,i)*deltaT*T[j/2][i]*c+ Src*deltaT*c;
     }
    // solve[0] += 0.5;
    }

    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
        //if (i == 10/2 || i == 31/2 || i == 100/2 || i == 316/2 || i == 1000/2 || i == 3162/2  || i == 10000/2) {
        if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000){
                if (j == 1 || j == 10 || j == 17 || j == 31 || j == 45 || j == 50 || j == 56 || j == 75 || j == 100 || j == 133 || j == 177)
                    printf("%f\t",E[j][i]);
               //     printf("%lf\t",pow(T[j][i],0.25));
            }
        }
            //if (i == 10/2 || i == 31/2 || i == 100/2 || i == 316/2 || i == 1000/2 || i == 3162/2  || i == 10000/2)
            if (i == 10 || i == 31 || i == 100 || i == 316 || i == 1000 || i == 3162  || i == 10000)
        printf("\n");
    }
    //printMatrix(E);
    sendToFile(E,T,deltaX,deltaT);
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

/*
  solves the equation Ax=v.
  x- initialy has the v vector, later it will have the x vector
  in our case Au(x,t+1) = u(x,t), that means x intiailly has u(x,t)
  U - the upper diagonal.
  L - the lower  diagonal
  mainD - the main Diagonal
*/
void solveTriagonal(double(*solve)[n],double L[n],double U[n],double mainD[n])
{
  int i,j;
  U[0] = U[0] / mainD[0];
  (*solve)[0] = (*solve)[0] / mainD[0];
  for (j = 0; j < n; j++) {
       // printf("%lf\t",mainD[j]);
        
    }
    // printf("\n");
  /* loop from 1 to X - 1 inclusive, performing the forward sweep */
  for (i = 1; i < n; i++) {
      const double m = 1.0f / (mainD[i] - L[i] * U[i - 1]);
 
      U[i] = U[i] * m;
           
      (*solve)[i] = ((*solve)[i] - (L[i] * ((*solve)[i - 1]))) * m;
     // printf("%lf\t", (*solve)[i]);
      
  }

  /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
  for (i = n - 2; i >=0 ; i--){
      (*solve)[i] -= U[i] * (*solve)[i + 1];
      // printf("BEFORE\tsolve: %15.15lf\tU: %15.15lf\ti is: %d\tcurrentsolve is:%15.15lf\n",(*solve)[i+1],U[i],i,(*solve)[i]);
        //(*solve)[i] -= (U[i] * (*solve)[i + 1]);

        // printf("AFTER\tsolve: %15.15lf\tU: %15.15lf\ti is: %d\tcurrentsolve is:%15.15lf\n",(*solve)[i+1],U[i],i,(*solve)[i]);
       }
  //  printf("\n\n");
}

/*
  builds the L,U,mainD arrays for the Solve triagonal
*/
void constructLUD(double (*L)[n],double (*U)[n],double (*mainD)[n],double mat[n][n]) {
  int i,j;
  for (i = 0; i < n-1; i++) {
    (*U)[i] = mat[i][i+1];
  }
    (*L)[0] = 0;
  for ( i = 1; i < n; i++) {
      (*L)[i] = mat[i][i-1];
  }
  for (i = 0; i < n; i++) {
    (*mainD)[i] = mat[i][i];
  }
}

//copies from solve to the matrix
void copyFromSolution(double*solve,double(*E)[X+1][N],double(*F)[X][N],int j) {
  int i = 0,k = 0;
  for ( i = 0; i < n; i+=2) {
      if (k == X) {
          break;
      }
    (*E)[k][j] = solve[i];
    (*F)[k][j] = solve[i+1];
    k++;
  }

  (*E)[X][j] = solve[n-1];
}

//copies from the matrix to the solve
void copyToSolution(double(*solve)[N],double mat[N][N],int j){
  int i;
  for ( i = 0; i < N; i++) {
    (*solve)[i] = mat[i][j];
  }
}

/**
 *
 * CONSIDERING A TRIAGONAL MATRIX ! ! !
*/
void matrixDotVector(double mat[][N], double (*vec)[N][N],int j) {
  int i,k;
  //vec[i][j+1] refers to E(x,t+1).
  for ( i = 0; i < N; i++) {
    for ( k = 0; k < N; k++) {
  //    (*vec)[i][j+1] += mat[i][k] *vec[k][j];
    }
  }
}

/*sets the triagonal matrix for E
lambda = deltaT, lambda2 = deltaX
*/
void setMatrixForE(double (*mat)[n][n],double deltaX,double deltaT,double sigmaT) {
  int i = 0;
  int j = 0;
  int k = 0;
  for (i = 1; i < n; i+=2) {
        (*mat)[i][i] = 1.0 + deltaT*getOpacity(k,currentTimeStep-1)*c;
        if (i != 0 ) {
          (*mat)[i][i-1] = -deltaT*c*c/(3.0*deltaX);
        }
       if (i != n-1) {
          (*mat)[i][i+1] = deltaT*c*c/(3.0*deltaX);
      }
      k++;
      printf("%lf\n",-deltaT*c*c/(3.0*deltaX));
  }
  k = 0;
  for (i = 0; i < n; i+=2) {
    (*mat)[i][i] = 1.0 + deltaT*getOpacity(k,currentTimeStep-1)*c ;
    if (i != 0) {
      (*mat)[i][i-1] = -deltaT/deltaX;
    }
    if (i != n-1) {
      (*mat)[i][i+1] = deltaT/deltaX;
    }
  }
  //(*mat)[0][0] += 0.5*lightspeed;
}


//prints the matrix
void printMatrix(double mat[][N]) {
  int i,j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%f\t", mat[i][j]);

    printf("\n");
  }
}

void printMatrix2(double mat[][n]) {
  int i,j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%f\t", mat[i][j]);

    printf("\n");
  }
}

void sendToFile(double E[N][N],double T[N][N],double deltaX,double deltaT) {
    int i = 0,j;
    FILE*fp;
    fp = fopen("../data/SuOlsonP1Data.txt","w");
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

double getOpacity(int space,int time1) {
  if (constOpacity) {
    return 1.0;
  } else {
      
        double v = T[space][time1];
        double t  = pow(v,0.25);
        double a = 1.0/(pow(t,3.0));
        return a;
        }
}

double getCv(int space, int time1) {
    if (constOpacity) {
        double a = getT(space,time1);
        return (alpha*pow(a,3));
    } else {
        return 4.0*arad;
    }
}

double getFinc(){
    return c*arad/4.0;
}

double getT(int space,int time1) {
    return pow(T[space][time1],0.25);
}