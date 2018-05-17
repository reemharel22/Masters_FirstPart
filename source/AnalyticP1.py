#this calculates the analytic solution of P1 Eq.
#THE PROBLEM IS, WHEN YOU TRY TO CALCULATE, ALTHOUGH HEAVISIDE GIVES 0, IT WILL CONTINUE THE WHOLE EQUATION !
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.special as special
import scipy
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.integrate import dblquad
#### THIS CALCULATES THE ANALYTIC SOLUTION ONLY FOR t = 1!
#### IF U WANT TO CHANGE IT, CHANGE PARAMETER t !!!
#### VERY EASY !
E = [];
I0 = [];
I1 = [];
f = open("../data/SuOlsonP1Data.txt","r")
y = []
#we only need the x values
with open("../data/SuOlsonP1Data.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers
        # (as strings)
        numbers_str = line.split();
        #convert numbers to floats
        numbers_float = [float(x) for x in numbers_str];
        y.append(numbers_float);
        break;

def integrand1(x,t):
    if (t > 10):
        a = integrate.quad(lambda y: np.exp(-math.sqrt(3)*math.fabs(x-y) ) *heaviside1(t- math.sqrt(3)*math.fabs(x-y))
        *heaviside1(10 - t + math.sqrt(3)*math.fabs(x-y) ),-0.5,0.5 );
    else:
        a = integrate.quad(lambda y: np.exp(-math.sqrt(3)*math.fabs(x-y) ) *heaviside1(t- math.sqrt(3)*math.fabs(x-y) ),-0.5,0.5 );
    b = (  a[0] *(math.sqrt(3) / 2)  )
    return b

def heaviside1(x):
    return 0.5 * (np.sign(x) + 1)

def f(y,T,x):
    if heaviside1(T- math.sqrt(3)*math.fabs(x-y)) != 0:
        a = math.sqrt( math.pow(T,2) - 3*math.pow((x-y),2)   )
        if a != 0:
            return (T * np.exp(-T)*heaviside1( T- math.sqrt(3)*math.fabs(x-y) ) * special.i1(a) ) / a;
    return 0

def sumXiC(n,k,x,t):
    sum = 0;
    if (t > 10) :
        for i in range(1,n-1):
            sum = sum + f(-0.5,t -10 + i*k,x) + f(0.5, t-10 +  i*k ,x)
    else :
        for i in range(1,n-1):
            sum = sum + f(-0.5,i*k,x) + f(0.5, i*k ,x)
    return 2*sum

def sumCXi(m,h,x,t):
    sum = 0
    if (t > 10) :
        for i in range(1,m-1):
            sum = sum + f(-0.5 + i*h ,t - 10,x) + f(-0.5 + i*h, t ,x)
    else :
        for i in range(1,m-1):
            sum = sum + f(-0.5 + i*h ,0,x) + f(-0.5 + i*h, t ,x)
    return 2*sum

def doubleSum(m,n,h,k,x,t):
    sum = 0
    if (t > 10):
        for i in range(1,n-1):
            for j in range(1,m-1):
                sum = sum + f(-0.5 + i*h ,t - 10 + j*k,x)
    else :
        for i in range(1,n-1):
            for j in range(1,m-1):
                sum = sum + f(-0.5 + i*h ,j*k,x)
    return 4*sum

def trap2d(x,t):
    m = 200
    n = 200;
    h = (1.0/m);
    k = t/n;
    if (t > 10):
        k = 10.0/n;
    s = 0.25*h*k * (  f(-0.5,0,x) + f(0.5,0,x) + f(-0.5,t,x) + f(0.5,t,x)
        + sumXiC(n,k,x,t) + sumCXi(m,h,x,t) + doubleSum(m,n,h,k,x,t)  );
    return s

tarr = [0.1, 0.316, 1.0, 3.16 , 10.0 , 31.6228, 100];

f2 = open('../data/P1AnalyticData.txt','w')
#sa = trap2d(0.13,t)*(math.sqrt(3) / 2)  + integrand1(0.12,t)
#print sa
for j in range(7):
    t = tarr[j];
    sa = "";
    for i in range(5000):
        x = (y[0][i])
        if (x == 0.5):
            x = x + 0.001;
        k = trap2d(x,t)*(math.sqrt(3) / 2)  + integrand1(x,t)
        if (k == 0):
            break;
        i = i + 1
        f2.write(str(k) + " ")
    f2.write("\n");
    print "done time - " + str(j);
f2.close
#this does the part where u calculate the dt part !!!!
