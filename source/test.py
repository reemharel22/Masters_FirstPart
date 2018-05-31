import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
#### THIS PLOTS THE BENCHMARK DIFFUSION/TRANSPORT AND MY SOLUTION !
#### TO CHANGE THE t ! (current t = 1) CHANGE IN LINE t0 TO THE WANTED TIME !
#### DONT FORGET TO CHANGE THE OTHER PYTHON FILE
#### WHAT ABOUT THE OTHER BENCHMARKS ?!?!
#### THE ONLY TIMES THAT ARE THERE ARE t = 3.16 t = 10 !!! ANY MORE THAN THAT YOU WILL
#### NEED TO ADD TO THE FILES : SuOlsonDiffusionData & SuOlsonTransportData
TT = 3;
t0 = 1
def SuOlsonMyNumericSolution(fname):
    k = 0;
    y = []
    i = 0;
    myList = [];
    f = open(fname,"r")
    with open(fname) as f:
        for line in f:  #Line is a string
            numbers_str = line.split()
            #convert numbers to floats
            numbers_float = [float(x) for x in numbers_str]
            myList.append(numbers_float)
            i = i + 1;
            if (i == int(t0*100) + 3):
                y.append(myList[ int(t0*100) +2])
                return y,myList[0]
    y.append(myList[ int(t0*100) +2])
    return y,myList[0]
def XXX(fname):
    k = 0;
    y = []
    i = 0;
    myList = [];
    f = open(fname,"r")
    with open(fname) as f:
        for line in f:  #Line is a string
            numbers_str = line.split()
            #convert numbers to floats
            numbers_float = [float(x) for x in numbers_str]
            myList.append(numbers_float)

    return myList
def P1MyNumericSolution(fname):
    k = 0;
    y = []
    myList = [];
    f = open(fname,"r")
    with open(fname) as f:
        for line in f:  #Line is a string
            numbers_str = line.split()
            #convert numbers to floats
            numbers_float = [float(x) for x in numbers_str]
            myList.append(numbers_float)
  #    for i in range(0,N):
    #    if (myList[1][i] == t0):
    y.append(myList[ int(t0*100) +2])
    return y,myList[0]
def GetFromFile(fname):
    f = open(fname,"r")
    y = []
    with open(fname) as f:
        for line in f:  #Line is a string
            #split the string on whitespace, return a list of numbers
            # (as strings)
            numbers_str = line.split()
            #convert numbers to floats
            numbers_float = [float(x) for x in numbers_str]
            y.append(numbers_float)
    return y;
#gets number of Y poitns, scatters them according to the transport X,Y
def ScatterdY(transportX,points,x1):
    y = []
    for j in range(0,11):
        for i in range(0,N):
            if (x1[i] == transportX[j]):
                y.append(points[i])
    return y;
#start of the main
N = 5001;
#MX = XXX("../data/MX_11.txt");
#MY = XXX("../data/MY_11.txt");
##LPX = XXX("../data/LPX_1.txt");
#LPY = XXX("../data/LPY_1.txt");
suOlsonP1Numerit = []
p1Analytic = []
#FLM,x1 = SuOlsonMyNumericSolution("../data/SuOlsonFluxLimitersDataMinerbo.txt");
#FLLP,x1 = SuOlsonMyNumericSolution("../data/SuOlsonFluxLimitersDataLP.txt");
#FLK,x1 = SuOlsonMyNumericSolution("../data/SuOlsonFluxLimitersDataKershaw.txt");
#p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#EFM,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo.txt");
#EFM2,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo2.txt");
suOlsonP1Numerit,x2 = P1MyNumericSolution("../data/SuOlsonP1Data.txt")
#EFLP,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorLP.txt");
#EFK,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorKershaw.txt");
transportX = [0.01,0.1,0.17,0.31,0.45,0.5,0.56,0.75,1.0,1.33,1.77 ]

#suOlsonTransp = [];
#suOlsonTransp = GetFromFile("../data/SuOlsonTransportData.txt")
#line3, = plt.plot(suOlsonTransp[0],suOlsonTransp[TT],'^k',label="Transport");
#line4 = plt.plot(MX,MY,'--r',label="Minerbo");
#line56 = plt.plot(LPX,LPY,'--g',label="LP");
line5, = plt.plot(x1[0:108],p1Analytic[TT-1][0:500],'--k',label="P1-Analytic");
line4, = plt.plot(x1[0:500],suOlsonP1Numerit[0][0:500],'k',label="P1-Numerit")
line4, = plt.plot(x1[0:500],EFM[0][0:500],'g',label="F-L Minerbo");
#line44, = plt.plot(x1[0:500],FLM[0][0:500],'--g',label="F-L Minerbo");
#line54, = plt.plot(x1[0:500],FLLP[0][0:500],'--r',label="F-L Levermore Pomraning");
#line64, = plt.plot(x1[0:500],FLK[0][0:500], '--b',label="F-L Kershaw");
#line44, = plt.plot(x1[0:500],EFM2[0][0:500],'k',label="F-L Minerbo - Avner");
line5, = plt.plot(x1[0:500],EFLP[0][0:500],'r',label="F-L Levermore Pomraning");
line6, = plt.plot(x1[0:500],EFK[0][0:500], 'b',label="F-L Kershaw");
#legend
x = [0.5,1,2,3,4]
y = [0,0.001,0.01,0.1,1,3]
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.3,4)
plt.ylim(1e-4,1.1)
#plt.axis([0.3,8,0.001,3]);
#ticks = [0.5,1,2,5,8];
#ticks2 = [0.001,0.01,0.1,1]
#labels = ['0.5','1','2','5','8']
#labels2 = ['0.001','0.01','0.1','1']
ticks = [0.3,0.5,1,1.4,1.6,3];
ticks2 = [0.001,0.01,0.1,1]
labels = ["0.3",'0.5','1','1.4','1.6','3']
labels2 = ['0.001','0.01','0.1','1']
plt.xticks(ticks,labels);
plt.yticks(ticks2,labels2)
plt.legend(prop={'size': 10})
plt.title("For t = " + str(t0))
plt.ylabel('Radiation Energy Density - E(x,t)');
plt.xlabel('x');
plt.show();
