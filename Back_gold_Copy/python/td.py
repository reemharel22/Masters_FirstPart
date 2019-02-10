import re
N = 1500
f = open('TD.data', 'r')
deltaT = [0] * 1500
Y = []
T = []
for i in range (0,1500):
    deltaT[i] = 0.002 * i
i = 0
with open('TD.data') as f:
    for line in f:
        a = re.split(', ', line)
        Y.append(float(a[1]))
        T.append(float(a[0]))
        i = i + 1

closest_value = 100 
closest_abs = 100
current_abs = 0
final_value = [0.0]*1500
for i in range (0,1500):
    closest_abs = 100
    for k in range(len(T)):
        closest_value = 0
        current_abs = abs(T[k] - deltaT[i])
        if closest_abs > current_abs:
            closest_abs = current_abs
            closest_value = Y[k]
            final_value[i] = closest_value
f = open('TD_1500_interpolate1.data', 'w+')
for i in range(len(final_value)):
    f.write(str(final_value[i]))
    f.write("\n")
f.close()