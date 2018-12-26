import os 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time


def get_th_t(f_name):
    f = open(f_name,"r")
    t = np.ndarray([759], np.float64)
    th = np.ndarray([759], np.float64)
    i = 0
    k = 0
    j = 0
    with open(f_name) as f:
        for line in f: 
            if i % 2 == 0:
                t[k] = float(line) * 1E-9
                k = k + 1
            else:
                th[j] = float(line) * 11605
                th[j] = 1.9 * 1160500
                j = j + 1
            i = i + 1
    return t, th

def sort_t_th(th, t, iters):
    for i in range(iters - 1):
        for j in range (iters - i - 1):
            if t[j] > t[j + 1]:
                tmp1 = t[j + 1]
                tmp2 = th[j + 1]
                t[j + 1] = t[j]
                th[j + 1] = th[j]
                t[j] = tmp1
                th[j] = tmp2
    return t, th

iterations = 759
# check if it's because of the sort ... maybe try with lesser points..
HeV = 1160500
alpha = 3.5
beta = 1.1
mu = 0.1
rho = 0.05
lambda_ = 0.75
g = 1.0/9175.0
f = 8.8E13
f = f / (pow(HeV, beta))
g = g / pow(HeV, alpha)
sigma_boltz = 5.6704E-5
four_alpha = 4.0 + alpha
epsilon = beta / four_alpha
print epsilon
c1 = 16 * g * sigma_boltz
c2 = 3.0 * f * pow(rho, 2.0 - mu + lambda_) * four_alpha
C  = c1 / c2 
T, TH = get_th_t("dataset1.csv")
T, TH = sort_t_th(TH,T,iterations)
x = np.ndarray([iterations])
th = TH[0]
H = pow(th, four_alpha)
H_integral = T[0] *  H
for i in range(1, 759):
    dt_curr = T[i]
    dt_prev = T[i - 1]
    dt = dt_curr - dt_prev
    if dt < 0:
        print "Negative dt..."
        print i
    th = TH[i]
    H = pow(th, four_alpha)
    H_integral = H * dt_curr
    epsilon_2 = 2.0 + epsilon
    p1 = epsilon_2 * C * pow(H, -epsilon)
    x[i] = p1 * H_integral / (1.0 - epsilon)

x = pow(x, 0.5)
x = x*1E1
ff = open("HS_analytic.data", "w+")
ff.writelines(str(T))
ff.writelines(str(x))
line441, = plt.plot(T / 1E-9,x,'k',label="P1AB")

plt.show()
