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
t0 = 3.16


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
#x1,t1,suOlsonWF, a1, flux = WaveFront("../data/Temp/Back_1500_WaveFront.txt")
#p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#EFM,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonEddingtonFactorMinerbo.txt");
#EFM2,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo2.txt");
#suOlsonDiffNumerit,x2 = SuOlsonMyNumericSolution("../data/Temp/diffsuindt1.txt")
#suOlsonDiffMUB,x1 =SuOlsonMyNumericSolution1("../data/Temp/Diff_TH5_DT02.txt")
#suOlsonP1Numerit,x1 = SuOlsonMyNumericSolution("../data/Temp/Diff_TH5_DT01.txt")
#p1Analytic,x1 = SuOlsonMyNumericSolution3("../data/Temp/Diff_TH5_DT001.txt")
plt.xlim(0, 2)
x2 , y1, y2, t1,t2 = dh.extract_data("../data/SuOlsonData.txt", t0)
line16, = plt.plot(x2, y1[0:2000],'g',label="weff")
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
plt.title("Su Olson for t = " + str(t1))
#plt.ylabel('Radation tempreture Density - T(x,t)');
#plt.xlabel('x');
plt.show()
