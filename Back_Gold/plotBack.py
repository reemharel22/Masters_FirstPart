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
NN = 1000
t0 = 0.01
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
y_half = []
y_one = []
y_two = []
y_three = []
y_third = []
y_2 = []
y_1 = []
tm_half = []
tm_one = []
tm_two = []
tm_three = []
tm_2 = []
tm_third = []

EV = 11605
x1,t1,suOlsonWF = dh.WaveFront("../data/Temp/Back_1500_WaveFront.txt")
#p1Analytic = GetFromFile("../data/P1AnalyticData.txt")
#EFM,x1 = SuOlsonMyNumericSolution("../data/Temp/SuOlsonEddingtonFactorMinerbo.txt");
#EFM2,x1 = SuOlsonMyNumericSolution("../data/SuOlsonEddingtonFactorMinerbo2.txt");
#suOlsonDiffNumerit,x2 = SuOlsonMyNumericSolution("../data/Temp/diffsuindt1.txt")
#suOlsonDiffMUB,x1 =SuOlsonMyNumericSolution1("../data/Temp/Diff_TH5_DT02.txt")
#suOlsonP1Numerit,x1 = SuOlsonMyNumericSolution("../data/Temp/Diff_TH5_DT01.txt")
#p1Analytic,x1 = SuOlsonMyNumericSolution3("../data/Temp/Diff_TH5_DT001.txt")
# x2 ,y_half,a1,t = dh.extract_data("../data/SuOlsonData.txt", 0.5)
# #x2 ,y_one,a1,t = dh.extract_data("../data/SuOlsonData.txt", 1)
# #x2 ,y_two,a1,t = dh.extract_data("../data/SuOlsonData.txt", 2)
# #x2 ,y_three,a1,t = dh.extract_data("../data/SuOlsonData.txt", t0)
# x2 ,y_third,a1,t = dh.extract_data("../data/SuOlsonData.txt", 0.25)
# x2 ,y_2,a1,t = dh.extract_data("../data/SuOlsonData.txt", 0.01)

# x1 ,tm_half,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", 0.5)
# #x1 ,tm_one,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", 1)
# #x2 ,tm_two,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", 2)
# #x2 ,tm_three,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", t0)
# x1 ,tm_third,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", 0.25)
# x1 ,tm_2,a1,t = dh.extract_data("../data/Temp/SuOlsonData.txt", 0.01)
x1, y1, y2, t1, t2 = dh.extract_data("../data/BackTemp.txt",t0)
x3, rad_1, rad_2, t1, t2 = dh.extract_data("../data/BackEnergy.txt",t0)
#plt.xlim(0,2)
for i in range(0,2000):
    rad_1[i] = rad_1[i] / EV
#    # x2[i] = x2[i] * 10
#     y_half[i] = y_half[i]/EV
#     #y_one[i] = y_one[i]/EV
#   #  y_two[i] = y_two[i]/EV
#    # y_three[i] = y_three[i]/EV
#     y_2[i] = y_2[i]/EV
#     y_third[i] = y_third[i]/EV
line16, = plt.plot(x1[0:2000], y1[0:2000])
line1, = plt.plot(x1[0:2000], rad_1[0:2000])
#line17, = plt.plot(t1,)
#line17, = plt.plot(x2[0:2000],y_one[0:2000],'g',label="TR t0 = 1")
#line18, = plt.plot(x2[0:2000],y_two[0:2000],label="TR t0 = 2")
#line19, = plt.plot(x2[0:2000],y_three[0:2000],label="TR t0 = 3")
#line29, = plt.plot(x2[0:2000],y_2[0:2000],'b',label="TR t0 = 0.25")
#line39, = plt.plot(x2[0:2000],y_third[0:2000],'k',label="TR t0 = 0.01")

#line16, = plt.plot(x1[0:2000],tm_half[0:2000],'--r',label="TM t0 = 0.5")#
#line17, = plt.plot(x1[0:2000],tm_one[0:2000],'--g',label="TM t0 = 1")
#line18, = plt.plot(x2[0:2000],tm_two[0:2000],label="TM t0 = 2")
#line19, = plt.plot(x2[0:2000],tm_three[0:2000],label="TM t0 = 3")
#line29, = plt.plot(x1[0:2000],tm_2[0:2000],'--b',label="TM t0 = 0.25")
#line39, = plt.plot(x1[0:2000],tm_third[0:2000],'--k',label="TM t0 = 0.01")
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
#line8, = plt.plot(x1[0:3000],suOlsonDiffMUB[0][0:3000],'--k',label="Disc Asym Diffusion")
#line4, = plt.plot(x5[0:3000],suOlsonDiffNumerit[0][0:3000],'b',label="Diffusion")
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
plt.title("For t = " + str(t1))
#plt.ylabel('Radation tempreture Density - T(x,t)');
#plt.xlabel('x');
plt.show()
