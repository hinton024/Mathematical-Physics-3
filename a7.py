import numpy as np
import matplotlib.pyplot as plt
# from IntegrationModule import *
from prettytable import PrettyTable
import numpy as np
from scipy.special.orthogonal import p_roots, l_roots, h_roots
from scipy import integrate

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


def MyLegQuadrature(f, a, b, n, m = 100):
    """
    Integrate `f` from `a` to `b` using  n-point composite Gauss Legendre Quadrature rule.

    Parameters
    ---------
    f : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    n : Integer 
        No. of points to integrate at.
    m : int, optional
        Number of subintervals for integration.
    
    Returns
    ---------
    leg : float
        n-point composite Gauss Legendre Quadrature approximation to the integral. 
    """
    h = (b - a) / m
    [leg_zer, w] = p_roots(n)
    leg = 0
    x_ = [a]
    for k in range(0, n):
        for i in range(1, m + 1):
            x_.append(a + h * i)
            leg += (h / 2) * w[k] * f(0.5 * h * leg_zer[k] + 0.5 * (x_[i] + x_[i - 1]))
    return leg

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

#First Representaition
def Representation1(x, epsilon, a = 0):
    return epsilon/(np.pi*((x-a)**2+epsilon**2))

#Second Representaition
def Representation2(x, epsilon, a = 0):
    return np.exp(-(x-a)**2/(2*epsilon))/(np.sqrt(2*np.pi*epsilon))


#Plots
x = np.linspace(-5, 5, 1000)


fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2)
for i in range(1, 6):
    ax1.plot(x, Representation1(x, 0.4/(2**i), a = 2), label = f"n = {i}")
    
ax1.set_xlim(1,3)
ax1.set_xlabel('x')
ax1.set_ylabel(r'$\phi_{\epsilon}(x)$')
ax1.set_title(r'$\phi_{\epsilon}(x)=\frac{\epsilon}{\pi((x-a)^2+\epsilon^2)}$ at a = 2')
ax1.legend()
ax1.grid()

for i in range(1, 6):
    ax2.plot(x, Representation1(x, 0.4/(2**i), a = -2), label = f"n = {i}")

ax2.set_xlim(-3,-1)
ax2.set_xlabel('x')
ax2.set_ylabel(r'$\phi_{\epsilon}(x)$')
ax2.set_title(r'$\phi_{\epsilon}(x)=\frac{\epsilon}{\pi((x-a)^2+\epsilon^2)}$ at a = -2')
ax2.legend()
ax2.grid()


fig, (ax3, ax4) = plt.subplots(nrows = 1, ncols = 2)
for i in range(1, 6):
    ax3.plot(x, Representation2(x, 0.4/(2**i), a = 2), label = f"n = {i}")

ax3.set_xlim(0,4)
ax3.set_xlabel('x')
ax3.set_ylabel(r'$\phi_{\epsilon}(x)$')
ax3.set_title(r'$\phi_{\epsilon}(x)=\frac{e^{\frac{-(x-a)^2}{2\epsilon}}}{\sqrt{2\pi\epsilon}}$ at a = 2')
ax3.legend()
ax3.grid()

for i in range(1, 6):
    ax4.plot(x, Representation2(x, 0.4/(2**i), a = -2), label = f"n = {i}")

ax4.set_xlim(-4,0)
ax4.set_xlabel('x')
ax4.set_ylabel(r'$\phi_{\epsilon}(x)$')
ax4.set_title(r'$\phi_{\epsilon}(x)=\frac{e^{\frac{-(x-a)^2}{2\epsilon}}}{\sqrt{2\pi\epsilon}}$ at a = -2')
ax4.legend()
ax4.grid()



plt.show()

#Integral 1
epsilon_list = [0.4/(2**1), 0.4/(2**2), 0.4/(2**3), 0.4/(2**4), 0.4/(2**5)]

i_val= []

legendre_11 = []
for i in range(1, 6):
    int = MyLegQuadrature(lambda x: Representation1(x, 0.4/2**i), -10, 10, 100)
    legendre_11.append(int)
    i_val.append(i)

hermite_11 = []
for i in range(1, 6):
    int = MyHermiteQuad(lambda x: Representation1(x, 0.4/2**i), 100)
    hermite_11.append(int)

simpson_11 = []
for i in range(1,6):
    int = My_Simp(lambda x: Representation1(x, 0.4/2**i), -10, 10, 10000)
    simpson_11.append(int)




legendre_12 = []
for i in range(1, 6):
    int = MyLegQuadrature(lambda x: Representation2(x, 0.4/2**i), -10, 10, 100)
    legendre_12.append(int)


hermite_12 = []
for i in range(1, 6):
    int = MyHermiteQuad(lambda x: Representation2(x, 0.4/2**i), 100)
    hermite_12.append(int)

simpson_12 = []
for i in range(1,6):
    int = My_Simp(lambda x: Representation2(x, 0.4/2**i), -10, 10, 10000)
    simpson_12.append(int)


table_11 = PrettyTable()
table_11.title = 'Integral I - First Representation'

table_11.field_names = ["n", "Legendre", "Hermite", "Simpson"]
for i in range(0,5):
    table_11.add_row([i_val[i], legendre_11[i], hermite_11[i], simpson_11[i]])

table_12 = PrettyTable()

table_12.field_names = ["n", "Legendre", "Hermite", "Simpson"]
table_12.title = 'Integral I - Second Representation'

for i in range(0,5):
    table_12.add_row([i_val[i], legendre_12[i], hermite_12[i], simpson_12[i]])



#Integral 2

legendre_21 = []
for i in range(1, 6):
    int = MyLegQuadrature(lambda x: Representation1(x, 0.4/2**i)*(x+1)**2, -10, 10, 100)
    legendre_21.append(int)

hermite_21 = []
for i in range(1, 6):
    int = MyHermiteQuad(lambda x: Representation1(x, 0.4/2**i)*(x+1)**2, 100)
    hermite_21.append(int)

simpson_21 = []
for i in range(1,6):
    int = My_Simp(lambda x: Representation1(x, 0.4/2**i)*(x+1)**2, -10, 10, 10000)
    simpson_21.append(int)

legendre_22 = []
for i in range(1, 6):
    int = MyLegQuadrature(lambda x: Representation2(x, 0.4/2**i)*(x+1)**2, -10, 10, 100)
    legendre_22.append(int)

hermite_22 = []
for i in range(1, 6):
    int = MyHermiteQuad(lambda x: Representation2(x, 0.4/2**i)*(x+1)**2, 100)
    hermite_22.append(int)

simpson_22 = []
for i in range(1,6):
    int = My_Simp(lambda x: Representation2(x, 0.4/2**i)*(x+1)**2, -10, 10, 10000)
    simpson_22.append(int)

table_21 = PrettyTable()
table_21.title = 'Integral II - First Representation'

table_21.field_names = ["n", "Legendre", "Hermite", "Simpson"]
for i in range(0,5):
    table_21.add_row([i_val[i], legendre_21[i], hermite_21[i], simpson_21[i]])

table_22 = PrettyTable()

table_22.field_names = ["n", "Legendre", "Hermite", "Simpson"]
table_22.title = 'Integral II - Second Representation'

for i in range(0,5):
    table_22.add_row([i_val[i], legendre_22[i], hermite_22[i], simpson_22[i]])



#Integral 3

legendre_31 = []
for i in range(1, 6):
    int = MyLegQuadrature(lambda x: Representation1(3*x+1, 0.4/2**i)*9*x**2, -10, 10, 100)
    legendre_31.append(int)

hermite_31 = []
for i in range(1, 6):
    int = MyHermiteQuad(lambda x: Representation1(3*x+1, 0.4/2**i)*9*x**2, 100)
    hermite_31.append(int)

simpson_31 = []
for i in range(1,6):
    int = My_Simp(lambda x: Representation1(3*x+1, 0.4/2**i)*9*x**2, -10, 10, 10000)
    simpson_31.append(int)


legendre_32 = []
for i in range(1, 6):
    int = MyLegQuadrature(lambda x: Representation2(3*x+1, 0.4/2**i)*9*x**2, -10, 10, 100)
    legendre_32.append(int)

hermite_32 = []
for i in range(1, 6):
    int = MyHermiteQuad(lambda x: Representation2(3*x+1, 0.4/2**i)*9*x**2, 100) 
    hermite_32.append(int)

simpson_32 = []
for i in range(1,6):
    int = My_Simp(lambda x: Representation2(3*x+1, 0.4/2**i)*9*x**2, -10, 10, 10000)
    simpson_32.append(int)

table_31 = PrettyTable()
table_31.title = 'Integral III - First Representation'

table_31.field_names = ["n", "Legendre", "Hermite", "Simpson"]
for i in range(0,5):
    table_31.add_row([i_val[i], legendre_31[i], hermite_31[i], simpson_31[i]])

table_32 = PrettyTable()

table_32.field_names = ["n", "Legendre", "Hermite", "Simpson"]
table_32.title = 'Integral III - Second Representation'

for i in range(0,5):
    table_32.add_row([i_val[i], legendre_32[i], hermite_32[i], simpson_32[i]])


print(table_11)
print(table_12)

print(table_21)
print(table_22)

print(table_31)
print(table_32)
 