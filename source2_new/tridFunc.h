#ifndef tridFunc_H_
#define tridFunc_H_

void solveTriagonal(int N,double(*solve)[N],double L[N],double U[N],double mainD[N]);
void constructLUD(int N,double (*L)[N],double (*U)[N-1],double (*mainD)[N],double mat[N][N]);
void printMatrix(int N,double mat[][N]);
void buildEFKershaw(int X,int N,double (*EF)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildEta(int X,int N,double (*EF)[X],double E[X][N],double F[X+1][N],int j);
void buildEFLP(int X,int N ,double (*EF)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildEFMinerbo(int X,int N,double (*EF)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildNoEF(int X,int N,double (*EF)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDKershaw(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDLP(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDMinerbo(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDDiff(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDDiffAsym(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDDiscDiffAsym(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double T[X][N],int j);
void setLambda(int X,int N,double E[][N],double T[][N],double (*D)[X],int j);
double getZ(double eta);
void BuildZ();

#endif
