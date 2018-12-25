import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import data_handle as dh

#### THIS PLOTS THE BENCHMARK DIFFUSION/TRANSPORT AND MY SOLUTION !
#### TO CHANGE THE t ! (current t = 1) CHANGE IN LINE t0 TO THE WANTED TIME !
#### DONT FORGET TO CHANGE THE OTHER PYTHON FILE
#### WHAT ABOUT THE OTHER BENCHMARKS ?!?!
#### THE ONLY TIMES THAT ARE THERE ARE t = 3.16 t = 10 !!! ANY MORE THAN THAT YOU WILL
#### NEED TO ADD TO THE FILES : SuOlsonDiffusionData & SuOlsonTransportData
TT = 4
t0 = 1
NN = 1000

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
x1,t1,suOlsonWF, a1, flux = dh.WaveFront("../data/Temp/Back_1500_WaveFront.txt")
#p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#EFM,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonEddingtonFactorMinerbo.txt");
#EFM2,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo2.txt");
#suOlsonDiffNumerit,x2 = SuOlsonMyNumericSolution("../data/Temp/diffsuindt1.txt")
#suOlsonDiffMUB,x1 =SuOlsonMyNumericSolution1("../data/Temp/Diff_TH5_DT02.txt")
#suOlsonP1Numerit,x1 = SuOlsonMyNumericSolution("../data/Temp/Diff_TH5_DT01.txt")
#p1Analytic,x1 = SuOlsonMyNumericSolution3("../data/Temp/Diff_TH5_DT001.txt")
x2 ,suOlsonDiffAsym,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", t0)
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
#plt.subplot(211)
plt.xlim(0,2)
line16, = plt.plot(x2[0:2000],suOlsonDiffAsym[0:2000],'g',label="weff");
#plt.subplot(212)
#line441, = plt.plot(t1[0][0:2000],suOlsonWF[0][0:2000],'k',label="P1AB")
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
plt.title("For t = " + str(t))
#plt.ylabel('Radation tempreture Density - T(x,t)');
#plt.xlabel('x');
plt.show()
