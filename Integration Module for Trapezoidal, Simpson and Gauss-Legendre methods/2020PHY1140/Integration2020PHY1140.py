from sympy import *
from scipy.special.orthogonal import p_roots
from scipy import integrate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import texttable as tt
tab = tt.Texttable()

x = symbols('x')
f = eval("lambda x:" + input("function to be integrated, f(x) = "))
#d = int(input("Enter the number of significant figures upto which results should be accurate = "))


def My_Trap(f_, a, b, n):
    y = []
    h = (b - a) / n
    for i in range(n + 1):
        y.append(f(a + i * h))  # y at limit points
    trp = h * (f(a) + f(b)) / 2
    for j in range(1, len(y) - 1):
        trp = trp + h * (y[j])
    return (trp)


def My_Simp(f_, a, b, n):
    h = (b - a) / (2 * n)
    simp = h * (f(a) + f(b)) / 3
    for i in range(1, 2 * n):
        if (i % 2 == 0):
            simp = simp + 2 * h * f(a + i * h) / 3
        elif (i % 2 == 1):
            simp = simp + 4 * h * f(a + i * h) / 3

    return (simp)


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
            print("integration using inbuilt function = ", integrate.quadrature(f, a, b))
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
            print("integration using inbuilt function = ", integrate.quadrature(f, a, b))
            break


def MyLegQuadrature(fs, a, b, n, m):  # Gauss legendre Quadrature ftion
    h = (b - a) / m
    e = []
    [leg_zer, w] = p_roots(n)
    leg_zer.tolist()
    w.tolist()
    sum_ = 0
    x_ = [a]
    s = []
    for k in range(0, n):
        for i in range(1, m + 1):
            x_.append(a + h * i)
            sum_ += (h / 2) * w[k] * f(0.5 * h * leg_zer[k] + 0.5 * (x_[i] + x_[i - 1]))
            s.append(sum_)
    return sum_


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


