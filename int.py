import numpy as np
from scipy.special.orthogonal import p_roots, l_roots, h_roots
from scipy import integrate


def My_Trap(f, a, b, n):
    """
    Integrate `f` from `a` to `b` using composite trapezoidal rule.
    
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
    trp : float
        Trapezoidal rule's approximation to the integral. 
    """
    y = []
    h = (b - a) / n
    for i in range(n + 1):
        y.append(f(a + i * h))
    trp = h * (f(a) + f(b)) / 2
    for j in range(1, len(y) - 1):
        trp = trp + h * (y[j])
    return trp


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


def MyLaguQuad(f, n):
    """
    Integrate `f` from `a` to `b` using  n-point  Gauss Laguerre Quadrature rule.

    Parameters
    ---------
    f : function
        A Python function or method to integrate.
    n : Integer 
        No. of points to integrate at.

    Returns
    ---------
    lagu : float
        n-point Gauss Laguerre Quadrature approximation to the integral. 
    """
    [lagu_zer, w] = l_roots(n)
    lagu = 0

    for i in range(1, n+1):
        lagu += f(lagu_zer[i-1])*w[i-1]
    
    return lagu


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



#Tolerance Functions

def MyTrap_tol(f_, a, b, n, d):
    i = 1
    err = []
    while i <= n:
        e = abs(f_("x", a, b, i) - f_("x", a, b, i + 1)) / f_("x", a, b, i + 1)
        err.append(e)
        if e > 0.5 * 10 ** -d:
            i = i + 1
            if i > n:
                print("Tolerance can't be reached for ", n, "Subintervals")
        elif e <= 0.5 * 10 ** -d:
            print("tolerance is reached in", i, "Intervals")
            print("integration using trapezoidal method = ", f_("x", a, b, i))
            print("integration using inbuilt function = ", integrate.quadrature(f_, a, b))
            break


def MySimp_tol(f_, a, b, n, d):
    i = 1
    err = []
    while i <= n:
        e = abs(f_("x", a, b, i) - f_("x", a, b, i + 1)) / f_("x", a, b, i + 1)
        err.append(e)
        if e > 0.5 * 10 ** -d:
            i = i + 1
            if i > n:
                print("Tolerance can't be reached for ", n, "intervals")
        elif e <= 0.5 * 10 ** -d:
            print("tolerance is reached in", i, "Intervals")
            print("integration using simpson method = ", f_("x", a, b, i))
            print("integration using inbuilt function = ", integrate.quadrature(f_, a, b))
            break

def MyLegQuadrature_tol(fs, a, b, n, m, d):
    i = 1
    err = []
    while i <= m:
        e = abs(MyLegQuadrature(fs, a, b, n, i) - MyLegQuadrature(fs, a, b, n, i + 1)) / MyLegQuadrature(fs, a, b, n,
                                                                                                         i + 1)
        err.append(e)
        if e > 0.5 * 10 ** -d:
            i = i * 2
            if i > m:
                print("Tolerance can't be reached for ", m, "Subintervals")
        elif e <= 0.5 * 10 ** -d:
            print("tolerance is reached in", i, "subintervals")
            print("integration and error using n point method(composite) = ", MyLegQuadrature(fs, a, b, n, i))
            # print("integration using inbuilt ftion = ", integrate.quadrature(f, a, b))
            break
            return MyLegQuadrature(fs, a, b, n, i)


def MyLegQuadrature_tol_1(fs, a, b, n, d):
    i = 1
    p = 0
    err = []
    while True:
        e = abs(MyLegQuadrature(fs, a, b, n, i) - MyLegQuadrature(fs, a, b, n, i + 1)) / MyLegQuadrature(fs, a, b, n,
                                                                                                         i + 1)
        err.append(e)
        if e > 0.5 * 10 ** -d:
            i = i * 2
        elif e <= 0.5 * 10 ** -d:
            p = MyLegQuadrature(fs, a, b, n, i)
            break
    return p, i, d


