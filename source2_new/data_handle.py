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

def extract_data(fname, t0):
    x = []
    y1 = []
    y2 = []
    dt = []
    i = 0
    f = open(fname,"r")
    data = f.readlines()
    x = convert_to_float(data[0].split())
    dt = convert_to_float(data[1].split())
    found = False
    while ( (i + 1) != len(dt) and not found):
        if (dt[i] == t0):
            y1 = convert_to_float(data[i + 2].split())
            y2 = y1
            t = dt[i]
            found = True
        if ( float(dt[i + 1]) > t0 and float(dt[i]) < t0):
            y1 = convert_to_float(data[i + 2].split())
            y2 = convert_to_float(data[i + 3].split())
            t = dt [i]

            found = True
        i = i + 1
    print ("Found i: ", i + 1)
    return x, y1, y2, t

def convert_to_float(arr):
    return [float(nm) for nm in arr]
