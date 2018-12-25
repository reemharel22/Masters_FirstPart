import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.integrate import dblquad

def create_TH(f_name, iter):
    fp1 = open(f_name,"r")
    TH = np.ndarray([iter])
    T = np.ndarray([iter])
    i = 0
    j = 0
    k = 0
    if iter == 759:
        with open(f_name) as f1:
            for line in f1:  #Line is a string
                if i % 2 == 0:
                    t = line
                    T[k] = float(t) * 1e-9
                    #T[k] = 0.01 * k * 1e-9
                    k = k + 1
                else:
                    y = line
                    TH[j] = float(y) * 11605
                    j = j + 1
                    i = i + 1
    i = 0
    if iter == 190:
        with open(f_name) as f1:
            for line in f1:  #Line is a string
                l = line.split(', ')
                l1 = l[1].split('\n')
                T[i] = float(l[0]) * 1e-9
                #TH[i] = float(l1[0]) * 11605
                TH[i] = 190 * 11605
                i = i + 1
    fp1.close()
    return T, TH

def create_file(T,TH) :
    f = open("dataset3.csv","w+")
    for i in range(0,190):
        f.writelines(str(T[i]/1e-9) + "\n")
        f.writelines(str(TH[i]/11605) + "\n")
    f.close()
        

HEV = 1160500
beta          = 1.1
alpha         = 3.5
lambda_       = 0.75
f             = 8.8e13 / ((HEV**beta))
g             = 1.0 / (9175.0 * (HEV**alpha))
rho = 50e-3
mu = 0.1

four_alpha = 4.0 + alpha
iterations = 190
epsilon = beta / four_alpha
sigma_boltzman = 5.670373e-5
C=16/(4+alpha)*g*sigma_boltzman/3/f/rho**(2-mu+lambda_)
TH = np.ndarray([iterations], np.float64)
T = np.ndarray([iterations], np.float64)
T, TH = create_TH("dataset2.csv", iterations)
x_f = np.ndarray([iterations], np.float64)
H2 = pow(TH[0], four_alpha)
eps_2 = 2 + epsilon
eps_1 = 1 - epsilon
create_file(T,TH)
H_integ = H2 * T[0]
k = 0


for i in range(1, iterations):    
    t = T[i]
    delta_t = T[i] - T[i - 1]
    H2 = pow(TH[i], four_alpha)
    H_integ = H_integ + H2 * (t - T[i - 1])
    if (t > 1 and T[i - 1] < 1):
        k = i - 1
    H1 = pow(TH[i], -four_alpha * epsilon)
    x_f[i -1] = (C * H1 * H_integ * (eps_2 / eps_1) )
    T1 = pow(TH[i], -four_alpha * epsilon)
    x_f[i-1] = pow(x_f[i - 1], 0.5) 

THY = np.ndarray([200], np.float64)
x = np.ndarray([200], np.float64)
plt.subplot(211)
line4 = plt.plot(T[0:iterations - 1]/1e-9, 10 * x_f[0:iterations - 1],'k',label="P1AB")
#from 0 to 0.2cm .... -> let's go for 200 points i.e 0.001cm

#delta_X = 0.001
#for i in range(0, 200):
#    x[i] = i * delta_X
#    x1 = 1.0 - (x[i] / x_f[k])
#    print (x1, x_f[k])
#    THY[i] = TH[k] * x1**(1/(four_alpha - beta))

#plt.subplot(212)
#line5 = plt.plot(x[0:iterations - 1] * 10, 100* THY[0:iterations - 1]/HEV,'k',label="P1AB")
plt.show()
