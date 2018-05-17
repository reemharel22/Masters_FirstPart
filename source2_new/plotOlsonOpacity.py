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
TT = 4;
t0 = 3;
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
suOlsonDiffNumerit = []
suOlsonP1Numerit = []
p1Analytic = []
suOlsonP1MUAB = []
suOlsonDiffMUB = []
suOlsonDiffAsym = []
#p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#EFM,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonEddingtonFactorMinerbo.txt");
#EFM2,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo2.txt");
#suOlsonDiffNumerit,x2 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonData.txt")
suOlsonP1Numerit,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonP1Data.txt")
#suOlsonDiffAsym,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionAsymptoticData.txt");
#suOlsonDiffMUB,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt");
#suOlsonP1MUAB,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonP1MUABData.txt")
#P1AB ,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonP1ABData.txt");
#transportX = [0.01,0.1,0.17,0.31,0.45,0.5,0.56,0.75,1.0,1.33,1.77 ]
#suOlsonDiffAsym,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionAsymptoticData.txt");
#suOlsonDiffMUB,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt");
#suOlsonTransp = [];
#suOlsonTransp = GetFromFile("../data/SuOlsonTransportData.txt")
#line3, = plt.plot(suOlsonTransp[0],suOlsonTransp[TT],'^g',label="Transport");
#line6, = plt.plot(x2[0:500],suOlsonP1MUAB[0][0:500],'k',label="P1 muAB");
line44, = plt.plot(x1[0:500],suOlsonP1Numerit[0][0:500],'b',label="P1")
#line7, = plt.plot(x1[0:500],suOlsonDiffAsym[0][0:500],'--b',label="Asym Diffusion")
#line8, = plt.plot(x2[0:500],suOlsonDiffMUB[0][0:500],'--k',label="Disc Asym Diffusion")
#line4, = plt.plot(x2[0:500],suOlsonDiffNumerit[0][0:500],'--b',label="Diffusion")
#line7, = plt.plot(x1[0:500],suOlsonDiffAsym[0][0:500],'--b',label="Asym Diffusion")
#line8, = plt.plot(x2[0:500],suOlsonDiffMUB[0][0:500],'--k',label="Disc Asym Diffusion")
#line5, = plt.plot(x1[0:500],P1AB[0][0:500],'b',label="P1-AB");
#legend
x = [0.5,1,2,3,4]
y = [0,0.001,0.01,0.1,1,3]
#plt.yscale('log')
#plt.xscale('log')
plt.ylim(0,1)
plt.xlim(0,2)
#plt.axis([0.3,8,0.001,3]);
#plt.xscale('log')
#ticks = [0.3,0.5,1,3,4];
#ticks2 = [0.001,0.01,0.1,1]
#labels = ['0.3','0.5','1','3','4']
#labels2 = ['0.001','0.01','0.1','1']
#plt.xticks(ticks,labels);
#plt.yticks(ticks2,labels2)
plt.legend(prop={'size': 10})
plt.title("For t = " + str(t0))
plt.ylabel('Radation tempreture Density - T(x,t)');
plt.xlabel('x');
plt.show();
