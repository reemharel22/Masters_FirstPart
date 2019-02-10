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
t0 = 1
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
x1, y1, y2, t1, t2 = dh.extract_data("../data/BackTemp.txt",t0)
x3, rad_1, rad_2, t1, t2 = dh.extract_data("../data/BackEnergy.txt",t0)
#plt.xlim(0,2)
for i in range(0, 1500):
    rad_1[i] = rad_1[i] / EV
    x1[i] = x1[i] * 10
line16, = plt.plot(x1[0:1500], y1[0:1500])
line1, = plt.plot(x1[0:1500], rad_1[0:1500])

plt.title("For t = " + str(t1))
#plt.ylabel('Radation tempreture Density - T(x,t)');
#plt.xlabel('x');
plt.show()
