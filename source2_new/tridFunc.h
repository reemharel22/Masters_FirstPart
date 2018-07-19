#ifndef tridFunc_H_
#define tridFunc_H_

void solveTriagonal(int N,float(*solve)[N],float L[N],float U[N],float mainD[N],float[N]);
void constructLUD(int N,float (*L)[N],float (*U)[N-1],float (*mainD)[N],float mat[N][N]);
void printMatrix(int N,float mat[][N]);
void buildEFKershaw(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildEta(int X,int N,float (*EF)[X],float E[X][N],float F[X+1][N],int j);
void buildEFLP(int X,int N ,float (*EF)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildEFMinerbo(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildNoEF(int X,int N,float (*EF)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDKershaw(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDLarsen(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDLP(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDMinerbo(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDDiff(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDDiffAsym(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float [X][N],int j);
void buildDDiscDiffAsym(int X,int N,float (*D)[X], float E[X][N],float F[X+1][N],float T[X][N],int j);
void setLambda(int X,int N,float E[][N],float T[][N],float (*D)[X],int j);
float getZ(float eta);
void BuildZ();

#endif
