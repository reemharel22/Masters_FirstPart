#ifndef SuOlsonAll_H_
#define SuOlsonAll_H_
float calculateWeff(int,int);
float calculateWeffNonAvg(int,int);
float getdx();
float getdt();
float calculateB(float);
float calculateMu2(float weff);
float getOpacity(int space,int time1);
float getFinc();
#endif
