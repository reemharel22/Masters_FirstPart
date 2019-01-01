#ifndef tridFunc_H_
#define tridFunc_H_

void solveTriagonal(int N,double(*solve)[N],double L[N],double U[N],double mainD[N]);
void constructLUD(int N,double (*L)[N],double (*U)[N-1],double (*mainD)[N],double mat[N][N]);
void printMatrix(int N,double mat[][N]);
void buildEFKershaw(int X,double (*EF)[X]);
void buildEta(int X,int N,double (*EF)[X],double E[X][N],double F[X+1][N],int j);
void buildEFLP(int X,double (*EF)[X]);
void buildEFMinerbo(int X,double (*EF)[X]);
void buildNoEF(int X,double (*EF)[X]);
void buildDKershaw(int X,double (*EF)[X]);
void buildDLarsen(int X,int N,double (*D)[X], double E[X][N],double F[X+1][N],double [X][N],int j);
void buildDLP(int X,double (*EF)[X]);
void buildDMinerbo(int X,double (*EF)[X]);
void buildDDiff(int X,double (*EF)[X]);
void buildDDiffAsym(int X,double (*EF)[X]);
void buildDDiscDiffAsym(int X,double (*EF)[X]);
void setLambda(int X,int N,double E[][N],double T[][N],double (*D)[X],int j);
double getZ(double eta);
void BuildZ();

#endif
