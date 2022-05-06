#Akarsh Shukla
#2020PHY1216
# Kabir sethi (2020PHY1097) and Brahmanand Mishra(2020PHY1184)


import matplotlib.pyplot as plt
from  MyIntegration  import  *
from matplotlib import  style
import pandas as pd
import numpy as np
from matplotlib import use
style.use("ggplot")
use("WebAgg")
def series(x, coeff,period):
    ser = []
    cosnx = lambda x,n: np.cos(n*np.pi*x/period)
    sinnx = lambda x,n: np.sin(n*np.pi*x/period)
    for i in x:
        temp = 0
        for j in range(len(coeff[0])):
            temp += cosnx(i,j)*coeff[0][j]
        for j in range(len(coeff[1])):
            temp += sinnx(i,j)*coeff[1][j]
        ser.append(temp)
    return ser

def fourier(f,  value,  N, d, period, method, fl = 1):
    a_n = []
    b_n = []
    meth_dict = {"simp":MySimp,"trap":MyTrap,"gauss":MyLegGauss}
    int_func = meth_dict[method]
    if value == 1 or value == 0:
        comb_func = lambda x,i: np.sin((i * np.pi * x) / period)
        for i in range(N):
            f1 = lambda x: comb_func(x,i) * f(x)
            if int_func == MyLegGauss:
                b_n.append(fl*(int_func(f = f1, a=0, b = period, n =12,m = 50)) / (period))
            else:
                a_n.append(method(fl*(int_func( a,  period, 1, f1,d)) / (period)))

    if value == -1 or value == 0:
        comb_func = lambda x,i: np.cos((i * np.pi * x) / period)
        for i in range(N):
            f1 = lambda x: comb_func(x,i) * f(x)
            if int_func == MyLegGauss:
                a_n.append(fl * (int_func(f=f1, a=0, b=period, n=12, m=50)) / (period))
            else:
                an.append(method(fl * (int_func(a, period, 1, f1, d)) / (period)))
            if i == 0:
                a_n[0] = a_n[0]/2
    return(a_n,b_n)



@np.vectorize
def sawtooth(x):
    if 0 <= x <=  np.pi:
        return x


@np.vectorize
def function1(x):
    if  type(x) != list:
        if -1 < x <= 0 or 1 < x < 2 or -2.5 <= x <= -2:
            return 0
        elif 0 < x <= 1 or -2 < x <= -1 or 2 < x <= 2.5:
            return 1
    else:
        temp = []
        for i in x:
            temp.append(function1(i))
        return temp

@np.vectorize
def function2(x):
    if  type(x) != list:
        if -1.5 <= x < -1 or 0.5 <= x < 1:
            return 0
        elif -2.5 <= x < -1.5 or -0.5 <= x < 0.5 or 1.5 <= x <= 2.5:
            return 1
        elif -1 <= x < -0.5 or 1 <= x < 1.5:
            return 0
    else:
        temp = []
        for i in x:
            temp.append(function2(i))
        return temp



@np.vectorize
def function3(x):
    if  type(x) != list:
        if -2.5 <= x < -2 or -1 <= x < 0 or 1 <= x < 2:
            return -0.5
        elif -2 <= x < -1 or 0 <= x < 1 or 2 <= x <= 2.5:
            return 0.5

    else:
        temp = []
        for i in x:
            temp.append(function2(i))
        return temp





list = [1, 2, 5, 10, 20]

print("FOR  X = -0.5")

x_check = [-0.5, 0, 0.5]
out_1 = []
for i in range(len(list)):
    t =  fourier(function1,0,  list[i],1e-10, 1, "gauss", 1)
    out_1.append((series(x_check, t, 1))[0])

comp_out_1 = function1(x_check)
rel = np.array(out_1) - comp_out_1[0]

dict = {"n values " : list, "series values" : out_1,  "relative error" : rel}
data = pd.DataFrame.from_dict(dict)
print(data)


print("FOR  X = 0")

out_2 = []

for i in range(len(list)):
    t =  fourier(function1,0,  list[i],1e-10, 1, "gauss", 1)
    out_2.append((series(x_check, t, 1))[1])

comp_out_1 = function1(x_check)
rel2 = np.array(out_2) - comp_out_1[1]

dict2 = {"n values " : list, "series values" : out_2,  "relative error" : rel2}
data2 = pd.DataFrame.from_dict(dict2)
print(data2)

print("FOR  X = 0.5")

out_3 = []

for i in range(len(list)):
    t =  fourier(function1,0,  list[i],1e-10, 1, "gauss", 1)
    out_3.append((series(x_check, t, 1))[2])

comp_out_1 = function1(x_check)
rel3 = np.array(out_3) - comp_out_1[2]

dict3 = {"n values " : list, "series values" : out_2,  "relative error" : rel3}
data3 = pd.DataFrame.from_dict(dict3)
print(data3)




x2 = np.linspace(-2.5, 2.5, 250)




value = []
for i in range(len(list)):
    t =  fourier(function1,0,  list[i],1e-10, 1, "gauss")
    value.append(series(x2, t, 1))
np.savetxt("Q1.txt", t)

figure, axis = plt.subplots(1,1)
plt.title("Q-1")
plt.plot(x2, value[0],label='n=1 terms',color="orange")
plt.plot(x2, value[1],label='n=2 terms',color="green")
plt.plot(x2, value[2],label='n=5 terms',color="yellow")
plt.plot(x2, value[3],label='n=10 terms',color="grey")
plt.plot(x2, value[4],label='n=20 terms',color="red")
plt.plot(x2, function1(x2))
plt.legend()
plt.savefig('q-1.pdf')



value2 = []
for i in range(len(list)):
    t =  fourier(function2,-1,  list[i],1e-10, 1, "gauss")
    value2.append(series(x2, t, 1))

figure, axis = plt.subplots(1,1)
plt.title("Q-2")
plt.plot(x2, value2[0],label='n=1 terms',color="orange")
plt.plot(x2, value2[1],label='n=2 terms',color="green")
plt.plot(x2, value2[2],label='n=5 terms',color="yellow")
plt.plot(x2, value2[3],label='n=10 terms',color="grey")
plt.plot(x2, value2[4],label='n=20 terms',color="red")
plt.plot(x2, function2(x2))
plt.legend()


value3 = []
for i in range(len(list)):
    t =  fourier(function3,0,  list[i],1e-10, 1, "gauss")
    value3.append(series(x2, t, 1))

figure, axis = plt.subplots(1,1)
plt.title("Q-3")
plt.plot(x2, value3[0],label='n=1 terms',color="orange")
plt.plot(x2, value3[1],label='n=2 terms',color="green")
plt.plot(x2, value3[2],label='n=5 terms',color="yellow")
plt.plot(x2, value3[3],label='n=10 terms',color="grey")
plt.plot(x2, value3[4],label='n=20 terms',color="red")
plt.plot(x2, function3(x2))
plt.legend()

x2 = np.linspace(-3*np.pi, 3*np.pi, 250)

value4 = []
for i in range(len(list)):
    t =  fourier(sawtooth,-1,  list[i],1e-10, np.pi, "gauss", fl = 2)
    value4.append(series(x2, t, np.pi))

figure, axis = plt.subplots(1,1)
plt.title("Cosine Series")
plt.plot(x2, value4[0],label='n=1 terms',color="orange")
plt.plot(x2, value4[1],label='n=2 terms',color="green")
plt.plot(x2, value4[2],label='n=5 terms',color="yellow")
plt.plot(x2, value4[3],label='n=10 terms',color="grey")
plt.plot(x2, value4[4],label='n=20 terms',color="red")
plt.plot(x2, sawtooth(x2))
plt.legend()

value5 = []
for i in range(len(list)):
    t =  fourier(sawtooth,1,  list[i],1e-10, np.pi, "gauss", fl = 2)
    value5.append(series(x2, t, np.pi))

figure, axis = plt.subplots(1,1)
plt.title("Sine Series")
plt.plot(x2, value5[0],label='n=1 terms',color="orange")
plt.plot(x2, value5[1],label='n=2 terms',color="green")
plt.plot(x2, value5[2],label='n=5 terms',color="yellow")
plt.plot(x2, value5[3],label='n=10 terms',color="grey")
plt.plot(x2, value5[4],label='n=20 terms',color="red")
plt.plot(x2, sawtooth(x2))
plt.legend()

print("FOR COSINE HALF SERIES")
x_check = [0, np.pi/2, np.pi]
out_1 = []
for i in range(len(list)):
    t =  fourier(sawtooth,-1,  list[i],1e-10, 1, "gauss", 1)
    out_1.append((series(x_check, t, 1))[0])

comp_out_1 = sawtooth(x_check)
rel = np.array(out_1) - comp_out_1[0]

dict = {"n values " : list, "series values" : out_1,  "relative error" : rel}
data = pd.DataFrame.from_dict(dict)
print(data)

out_2 = []

for i in range(len(list)):
    t =  fourier(sawtooth,-1,  list[i],1e-10, 1, "gauss", 1)
    out_2.append((series(x_check, t, 1))[1])

comp_out_1 = sawtooth(x_check)
rel2 = np.array(out_2) - comp_out_1[1]

dict2 = {"n values " : list, "series values" : out_2,  "relative error" : rel2}
data2 = pd.DataFrame.from_dict(dict2)
print(data2)



out_3 = []

for i in range(len(list)):
    t =  fourier(sawtooth,-1,  list[i],1e-10, 1, "gauss", 1)
    out_3.append((series(x_check, t, 1))[2])

comp_out_1 = sawtooth(x_check)
rel3 = np.array(out_3) - comp_out_1[2]

dict3 = {"n values " : list, "series values" : out_2,  "relative error" : rel3}
data3 = pd.DataFrame.from_dict(dict3)
print(data3)



print("FOR SINE HALF SERIES")
x_check = [0, np.pi/2, np.pi]
out_1 = []
for i in range(len(list)):
    t =  fourier(sawtooth,1,  list[i],1e-10, 1, "gauss", 2)
    out_1.append((series(x_check, t, 1))[0])

comp_out_1 = sawtooth(x_check)
rel = np.array(out_1) - comp_out_1[0]

dict = {"n values " : list, "series values" : out_1,  "relative error" : rel}
data = pd.DataFrame.from_dict(dict)
print(data)

out_2 = []

for i in range(len(list)):
    t =  fourier(sawtooth,1,  list[i],1e-10, 1, "gauss", 2)
    out_2.append((series(x_check, t, 1))[1])

comp_out_1 = sawtooth(x_check)
rel2 = np.array(out_2) - comp_out_1[1]

dict2 = {"n values " : list, "series values" : out_2,  "relative error" : rel2}
data2 = pd.DataFrame.from_dict(dict2)
print(data2)



out_3 = []

for i in range(len(list)):
    t =  fourier(sawtooth,1,  list[i],1e-10, 1, "gauss",  fl = 2)
    out_3.append((series(x_check, t, 1))[2])

comp_out_1 = sawtooth(x_check)
rel3 = np.array(out_3) - comp_out_1[2]

dict3 = {"n values " : list, "series values" : out_2,  "relative error" : rel3}
data3 = pd.DataFrame.from_dict(dict3)
print(data3)
plt.show()






