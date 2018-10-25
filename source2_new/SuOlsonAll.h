#ifndef SuOlsonAll_H_
#define SuOlsonAll_H_
double calculateWeff(int,int);
double calculateWeffNonAvg(int,int);
double getdx();
double getdt();
double calculateB(double);
double calculateMu2(double weff);
double getOpacity(int space,int time1);
double getFinc();
#endif
