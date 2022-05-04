
import math
import numpy as np
from scipy.special import eval_legendre
import matplotlib.pyplot as plt

def MySimp(a, b, n, f, Pnx, tbi,p):
  h = (b-a)/n 
  x_arr = []
  for i in range(0,n+1):
    x_ele = a + i*h 
    x_arr.append(x_ele)
  y_arr = []
  for i in range(0, n+1): 
    y_ele = tbi(f,Pnx,p,x_arr[i]) 
    y_arr.append(y_ele)
  sum = 0
  for i in range(0, n+1):
    if i == 0 or i == n:
      sum = sum + y_arr[i]
    elif i % 2 == 0:
      sum = sum + 2*y_arr[i]
    elif i % 2 == 1:
      sum = sum + 4*y_arr[i]
  sum = sum*(h/3)
  return sum 

def f2(x):
  return math.sin(x)*math.cos(x)

def f1(x):
  return 2*x**4 + 3*x + 2

def Pnx(n,x):
  return eval_legendre(n,x)

def tbi(f,Pnx,n,x):
  return f(x)*Pnx(n,x)

def expand(f,terms):
  cn = []
  for i in range(0,terms):
    coeff = ((2*i+1)/2)*MySimp(-1, 1, 100, f, Pnx, tbi, i)
    cn.append(coeff)
  return cn

terms1 = 5
Pnxs = ["P0(x)","P1(x)","P2(x)","P3(x)","P4(x)","P5(x)","P6(x)","P7(x)","P8(x)","P9(x)","P10(x)"]
arr1 = expand(f1,terms1)
print("Non-zero terms in legendre series expansion of 2*x**4 + 3*x + 2 : ")
for i in range(0,terms1):
  if round(arr1[i],2)!=0:
    if i<=3:
      print(f"{round(arr1[i],2)}*{Pnxs[i]}", end=" + ")
    else:
      print(f"{round(arr1[i],2)}*{Pnxs[i]}")
print("Coefficients of Pn(x) have been rounded off to 2 decimal places")

print("\nFirst 10 terms in legendre series expansion of sin(x)*cos(x): ")
terms2 = 10
arr2 = expand(f2,terms2)
for i in range(0,terms2):
  if round(arr2[i],4)<(-10e-15):
    if i<=8:
      print(f"({round(arr2[i],4)}*{Pnxs[i]})", end=" + ")
    else:
      print(f"({round(arr2[i],4)}*{Pnxs[i]})")
  else:
    if i<=8:
      print(f"{round(arr2[i],4)}*{Pnxs[i]}", end=" + ")
    else:
      print(f"{round(arr2[i],4)}*{Pnxs[i]}")
print("Coefficients of Pn(x) have been rounded off to 4 decimal places")
print()

fig1, (ax11, ax12) = plt.subplots(nrows=1, ncols=2, figsize=(15,5))
fig1.suptitle("Graphical Analysis of given functions")

ts = [1,2,3,4,5]
y_values1_1 = []
y_values1_2 = []
y_values1_3 = []
y_values1_4 = []
y_values1_5 = []
y_values1 = [y_values1_1,y_values1_2,y_values1_3,y_values1_4,y_values1_5]
err_arr = []

xs = np.linspace(-2,2,101)
py_val1_arr = []
for i in range(0,101):
  py_val1 = 2*xs[i]**4 + 3*xs[i] + 2
  py_val1_arr.append(py_val1)

for k in range(0,5):
  cn = expand(f1,ts[k])
  for i in range(0,101):
    expansion_value_at_xi = 0
    for j in range(0,ts[k]):
      expansion_value_at_xi = expansion_value_at_xi + cn[j]*Pnx(j,xs[i])
    y_values1[k].append(expansion_value_at_xi)

for i in range(0,5):
  ax11.plot(xs,y_values1[i], label=f"for n={i+1}")
ax11.plot(xs,py_val1_arr,label=f"inbuilt value", color='k')
ax11.grid()
ax11.set_xlabel("x")
ax11.set_ylabel(r'$2x^4 + 3x +2$')
ax11.set_title(r'Graphical analysis of $2x^4 + 3x +2$')
ax11.legend()

print()

ts = [2,4,6,8,10]
y_values2_1 = []
y_values2_2 = []
y_values2_3 = []
y_values2_4 = []
y_values2_5 = []
y_values2 = [y_values2_1,y_values2_2,y_values2_3,y_values2_4,y_values2_5]
err_arr = []

xs = np.linspace(-math.pi,math.pi,101)
py_val2_arr = []
for i in range(0,101):
  py_val2 = math.sin(xs[i])*math.cos(xs[i])
  py_val2_arr.append(py_val2)

for k in range(0,5):
  cn = expand(f2,ts[k])
  for i in range(0,101):
    expansion_value_at_xi = 0
    for j in range(0,ts[k]):
      expansion_value_at_xi = expansion_value_at_xi + cn[j]*Pnx(j,xs[i])
    y_values2[k].append(expansion_value_at_xi)

for i in range(0,5):
  ax12.plot(xs,y_values2[i], label=f"for n={2*(i+1)}")
ax12.plot(xs,py_val2_arr,label=f"inbuilt value", color='k')
ax12.grid()
ax12.set_xlabel("x")
ax12.set_ylabel(r'$sin(x)\cdot cos(x)$')
ax12.set_ylim(-1,1)
ax12.set_title(r'Graphical analysis of $sin(x)\cdot cos(x)$')
ax12.legend()
fig1.savefig("graph.pdf")