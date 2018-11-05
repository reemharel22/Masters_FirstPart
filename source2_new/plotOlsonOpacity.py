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
t0 = 1
def WaveFront(fname):
    f = open(fname,"r")

    x1 = []
    
    tt = f.readline().split()
    x1.append([float(x) for x in tt])


    t1 = []
    tt = f.readline().split()
    t1.append([float(x) for x in tt])


    data = f.readline()
    n_str = data.split()
    List = []
    List.append([float(x) for x in n_str])

    data = f.readline()
    n_str = data.split()
    a1 = []
    a1.append([float(x) for x in n_str])
    

    data = f.readline()
    n_str = data.split()
    flux = []
    flux.append([float(x) for x in n_str])
    
    return x1, t1, List, a1 ,flux

def SuOlsonMyNumericSolution(fname):
    k = 0
    y = []
    i = 0
    myList = []
    f = open(fname,"r")
    with open(fname) as f:
        for line in f:  #Line is a string
            numbers_str = line.split()
            #convert numbers to floats
            numbers_float = [float(x) for x in numbers_str]
            myList.append(numbers_float)
            i = i + 1;
            if (i == int(t0*1000) + 3):
                y.append(myList[ int(t0*1000) +2])
                return y, myList[0]
    y.append(myList[ int(t0*1000) +2])
    return y,myList[0]

def SuOlsonMyNumericSolution1(fname):
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
                if (i == int(t0*50) + 3):
                    y.append(myList[ int(t0*50) +2])
                    return y,myList[0]
        y.append(myList[ int(t0*50) +2])
        return y,myList[0]
def SuOlsonMyNumericSolution2(fname):
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
                if (i == int(t0*200) + 3):
                    y.append(myList[ int(t0*200) +2])
                    return y,myList[0]
        y.append(myList[ int(t0*200) +2])
        return y,myList[0]
def SuOlsonMyNumericSolution3(fname):
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
                if (i == int(t0*1000) + 3):
                    y.append(myList[ int(t0*1000) +2])
                    return y,myList[0]
        y.append(myList[ int(t0*1000) +2])
        return y,myList[0]
def SuOlsonMyNumericSolution4(fname):
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
                if (i == int(t0*400) + 3):
                    y.append(myList[ int(t0*400) +2])
                    return y,myList[0]
        y.append(myList[ int(t0*400) +2])
        return y,myList[0]
def AAA(fname):
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
def SuOlsonMyNumericSolution5(fname):
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
                if (i == int(t0*20) + 3):
                    y.append(myList[ int(t0*20) +2])
                    return y,myList[0]
        y.append(myList[ int(t0*20) +2])
        return y,myList[0]
def SuOlsonMyNumericSolution6(fname):
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
                if (i == int(t0*10) + 3):
                    y.append(myList[ int(t0*10) +2])
                    return y,myList[0]
        y.append(myList[ int(t0*10) +2])
        return y,myList[0]
#gets number of Y poitns, scatters them according to the transport X,Y
def ScatterdY(transportX,points,x1):
    y = []
    for j in range(0,11):
        for i in range(0,N):
            if (x1[i] == transportX[j]):
                y.append(points[i])
    return y;
#start of the main
N = 5001
suOlsonDiffNumerit = []
suOlsonP1Numerit = []
p1An1alytic = []
suOlsonP1MUAB = []
suOlsonDiffMUB = []
suOlsonDiffAsym = []
suOlsonWF = []
flux = []
a1 = []
x1,t1,suOlsonWF, a1, flux = WaveFront("../data/Temp/Back_1500_WaveFront.txt")
#p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#EFM,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonEddingtonFactorMinerbo.txt");
#EFM2,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo2.txt");
#suOlsonDiffNumerit,x2 = SuOlsonMyNumericSolution("../data/Temp/diffsuindt1.txt")
#suOlsonDiffMUB,x1 =SuOlsonMyNumericSolution1("../data/Temp/Diff_TH5_DT02.txt")
#suOlsonP1Numerit,x1 = SuOlsonMyNumericSolution("../data/Temp/Diff_TH5_DT01.txt")
#p1Analytic,x1 = SuOlsonMyNumericSolution3("../data/Temp/Diff_TH5_DT001.txt")
suOlsonDiffAsym,x2 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonData.txt")
#t0 = 20;
#f = open("dataset.csv","r")
#f1 = open("dataset1.csv","w+")
#data1 = f.readlines()
#f.close()
#f = open("dataset.csv","r")
#X = []
#for i in range(len(data1)):
#    data = f.readline()
#    x = float(data.split(',')[0])
#    y = float(data.split(',')[1])
#    X.append(str(x) + "\n" + str(y) + "\n")
#f1.writelines(X)
#f1.close()
#f.close()
#p1Analytic,x2 = SuOlsonMyNumericSolution("../data/Temp/001.txt")
#t0=30line441, = plt.plot(x1[0:3000],suOlsonDiffMUB[0][0:3000],'k',label="Diffusion dt=0.002")
#suOlsonP1Numerit,x1 = SuOlsonMyNumericSolution1("../data/Temp/002.txt")
#suOlsonP1Numerit,x4 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonP1Data.txt")
#suOlsonDiffMUB,x2 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt")
#suOlsonDiffAsym,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionAsymptoticData.txt");
#suOlsonDiffMUB,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt");
#suOlsonDiffAsym,x1 = SuOlsonMyNumericSolution("../data/SuOlsonData.txt")
#suOlsonP1AB,x2 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonP1ABData.txt")
#wef,x1 = SuOlsonMyNumericSolution("../data/Temp/weff.txt")
#P1 ,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonP1Data.txt")
#suOlsonDiffNumerit ,x5 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonData.txt");

#transportX = [0.01,0.1,0.17,0.31,0.45,0.5,0.56,0.75,1.0,1.33,1.77 ]
#suOlsonDiffAsym,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionAsymptoticData.txt");
#suOlsonDiffMUB,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonDiffusionDiscAsymptoticData.txt");
#suOlsonTransp = [];
#suOlsonTransp = GetFromFile("../data/SuOlsonTransportData.txt")
#line3, = plt.plot(suOlsonTransp[0],suOlsonTransp[TT],'^g',label="Transport");
#line6, = plt.plot(x3[0:3000],suOlsonP1MUAB[0][0:3000],'g',label="P1 MUAB")
#line6, = plt.plot(x2[0:3000],suOlsonP1AB[0][0:3000],'r',label="P1 AB")
plt.subplot(211)
line16, = plt.plot(x2[0:2000],suOlsonDiffAsym[0][0:2000],'g',label="weff");
plt.subplot(212)
line441, = plt.plot(t1[0][0:2000],suOlsonWF[0][0:2000],'k',label="P1AB")
#line441, = plt.plot(x5[0:1500],suOlsonDiffNumerit[0][0:1500],'k',label="P1AB")

#line44, = plt.plot(x5[0:1500], a1[0][0:1500],'r',label="Diffusion")
#line44, = plt.plot(t1[0][0:1500], flux[0][0:1500],'b',label="Diffusion")
#line4, = plt.plot(x1[0:1500], flux[0][0:1500],'g',label="Diffusion")
#a2 = [None] * 1500
#for i in range(1500):
#    a2[i] = a1[0][i] + flux[0][i]

#line4, = plt.plot(x1[0:1500], a2,'b',label="Diffusion")

#line7, = plt.plot(x1[0:500],suOlsonDiffAsym[0][0:500],'--b',label="Asym Diffusion")
#line8, = plt.plot(x1[0:3000],suOlsonDiffMUB[0][0:3000],'r',label="Disc Asym Diffusion")

#line411, = plt.plot(x1[0:3000],suOlsonDiffAsym[0][0:3000],'r',label="Diffusion dt = 0.005")
#line7, = plt.plot(x1[0:3000],p1Analytic[0][0:3000],'b',label="Diffusion dt = 0.001")
#line4, = plt.plot(x5[0:3000],suOlsonDiffNumerit[0][0:3000],'b',label="Diffusion")
#line8, = plt.plot(x1[0:3000],suOlsonDiffMUB[0][0:3000],'--k',label="Disc Asym Diffusion")
#line5, = plt.plot(x1[0:500],P1AB[0][0:500],'b',label="P1-AB");
#legend
#x = [0.5,1,2,3,4]
#y = [0,0.001,0.01,0.1,1,3]
#plt.yscale('log')
#plt.xscale('log')
#plt.xlim(0,20)
#plt.ylim(0,1)
#plt.axis([0.3,8,0.001,3]);
#plt.xscale('log')
#ticks = [0.3,0.5,1,3,4];
#ticks2 = [0.001,0.01,0.1,1]
#labels = ['0.3','0.5','1','3','4']
#labels2 = ['0.001','0.01','0.1','1']
#plt.xticks(ticks,labels);
#plt.yticks(ticks2,labels2)
#plt.legend(prop={'size': 10})
#plt.title("For t = " + str(t0))
#plt.ylabel('Radation tempreture Density - T(x,t)');
#plt.xlabel('x');
plt.show()
