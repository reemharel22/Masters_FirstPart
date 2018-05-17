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
t0 = 1;
def SuOlsonMyNumericSolution(fname):
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
    y.append(myList[ int(t0*50) +2])
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
mySuOlson = [];
p1Analytic = [];
scatterdMyDiffusion = [];
scatterdP1Analytic = []
scatterdP1Numerit = []
suOlsonDiff = []
x1 = []
x2 = []
transportX = [0.01,0.1,0.17,0.31,0.45,0.5,0.56,0.75,1.0,1.33,1.77 ]
suOlsonP1Numerit = []
suOlsonTransp = [];
suOlsonP1Numerit,x1 = P1MyNumericSolution("../data/SuOlsonP1Data.txt")
suOlsonTransp = GetFromFile("../data/SuOlsonTransportData.txt")
suOlsonDiff = GetFromFile("../data/SuOlsonDiffusionData.txt")
mySuOlson,x2 = SuOlsonMyNumericSolution("../data/SuOlsonData.txt")
p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#THIS IS TYHE WRONG DATA FROM THE p1 NUMERIT !!!!! NEED TO GO TO DIFFERENT LINE
#the 2 first items in myList will always be the x and t respectfuly
line1, = plt.plot(x1[0:500],mySuOlson[0][0:500],label="My Su-Olson")
line2, = plt.plot(suOlsonDiff[0],suOlsonDiff[7],'o',label="Su-Olson Diffusion")
line3, = plt.plot(suOlsonTransp[0],suOlsonTransp[7],'g^',label="Transport");
line4, = plt.plot(x2[0:500],suOlsonP1Numerit[0][0:500],'r',label="P1-Numerit")
line5, = plt.plot(x1[0:69],p1Analytic[2][0:69],'--b',label="P1-Analytic");
#legend

plt.legend(prop={'size': 14})
plt.xlim(0,5)
plt.ylim(0,0.4)
plt.title("For t = " + str(t0))
plt.show();
