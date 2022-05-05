#Ankur Kumar 2020PHY1113
#Preetpal Singh 2020PHY1140


# from IntegrationModule import *
from scipy.integrate import quad
import numpy as np
import pandas as pd
from prettytable import PrettyTable
from scipy.special.orthogonal import h_roots
def My_Simp(f, a, b, n):
    """
    Integrate `f` from `a` to `b` using composite simpson rule.

    Parameters
    ---------
    f : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    n : int
        Number of subintervals for integration.
    
    Returns
    ---------
    simp : float
        Simpson rule's approximation to the integral. 
    """
    h = (b - a) / (2 * n)
    simp = h * (f(a) + f(b)) / 3
    for i in range(1, 2 * n):
        if (i % 2 == 0):
            simp = simp + 2 * h * f(a + i * h) / 3
        elif (i % 2 == 1):
            simp = simp + 4 * h * f(a + i * h) / 3

    return simp

def MyHermiteQuad(f, n):
    """
    Integrate `f` from `a` to `b` using  n-point  Gauss Hermite Quadrature rule.

    Parameters
    ---------
    f : function
        A Python function or method to integrate.
    n : Integer 
        No. of points to integrate at.

    Returns
    ---------
    herm : float
        n-point Gauss Hermite Quadrature approximation to the integral. 
    """
    [herm_zer, w] = h_roots(n)
    herm = 0

    for i in range(1, n+1):
        herm += f(herm_zer[i-1])*w[i-1]
    
    return herm




ver_f = eval("lambda x:" + input("function to be integrated by Gauss Hermite , f(x) = "))

x = PrettyTable()

x.field_names = ["n-point", "Inbuilt Function", "My Function"]

x.add_row(["2", quad(lambda x: np.exp(-(x**2))*ver_f(x), -np.inf, np.inf)[0], MyHermiteQuad(ver_f, 2)])
x.add_row(["4", quad(lambda x: np.exp(-(x**2))*ver_f(x), -np.inf, np.inf)[0], MyHermiteQuad(ver_f, 4)])

print(x)

def Int_1(x):
    return 1/(1+x**2)

def Int_2(x):
    return np.exp(x**2)/(1+x**2)

I1 = []
I2 = []
n = []
for i in range(1, 8):
    val_1 = MyHermiteQuad(Int_1, 2**i)
    I1.append(val_1)


for i in range(1, 8):
    val_2 = MyHermiteQuad(Int_2, 2**i)
    I2.append(val_2)
    n.append(2**i)

data = {

  'n': n,
  'I_1': I1,
  'I_2': I2

}

df = pd.DataFrame(data)
numpy_array = df.to_numpy()
np.savetxt("quad-herm-1113.out.txt", numpy_array, fmt = "%f")
print(df)


#Comparsion With Simpson Method

func_1 = lambda x: np.exp(-(x**2))/(1+x**2)
func_2 = lambda x: 1/(1+x**2)

simp1 = []
simp2 = []
upperLimit = []
lowerLimit = []

for i in range(0,5):

    temp1 = My_Simp(func_1, -10**i, 10**i, 100000)
    simp1.append(temp1)
    
    temp2 = My_Simp(func_2, -10**i, 10**i, 100000)
    simp2.append(temp2)

    temp3 = 10**i
    upperLimit.append(temp3)

    temp4 = -10**i
    lowerLimit.append(temp4)


data1 = {
    'Lower Limit': lowerLimit,
    'Upper Limit': upperLimit,
  'I_1(Simpson)': simp1,
  'I_2(Simpson)': simp2

}

df = pd.DataFrame(data1)
print(df)
