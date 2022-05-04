import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from MyIntegration import MyLegQuadrature_tol_1, MyLegQuadrature, MyTrap_tol, MySimp_tol
import pandas as pd


x = sym.symbols('x')


def product(f1, f2):
    return lambda x: f1(x) * f2(x)


def sumfn(f1, f2):
    return lambda x: f1(x) + f2(x)


def FourierCoeff(method, d, f, dis, L, N, var):

    def coef_a0():

        a0 = 0
        low = dis[0]

        if method == 'quad':  # if method is quad

            for i in range(len(f)):
                # defining high index
                high = dis[i + 1]

                a0 += MyLegQuadrature_tol_1(f[i], low, high, n=5, d=d)
                low = high

        elif method == 'trap':

            for i in range(len(f)):
                # defining high index
                high = dis[i + 1]

                a0 += MyTrap_tol(f[i], low, high, n=20, d=d)
                low = high

        elif method == 'simp':

            for i in range(len(f)):
                # defining high index
                high = dis[i + 1]

                a0 += MySimp_tol(f[i], low, high, n=20, d=d)
                low = high

        return a0 / L

    def coef_an():
        an = []
        if method == 'quad':

            for i in range(1, N + 1):
                # i th coefficient of fourier series
                ai = 0

                f_for_a = lambda x: np.cos(i * np.pi * x / L)

                low = dis[0]

                for j in range(len(f)):
                    high = dis[j + 1]
                    ai += MyLegQuadrature_tol_1(product(f[j], f_for_a), low, high, n=5, d=d)

                    low = high

                an.append(ai)

        elif method == 'trap':

            for i in range(1, N + 1):
                # i th coefficient of fourier series
                ai = 0

                f_for_a = lambda x: np.cos(i * np.pi * x / L)

                low = dis[0]

                for j in range(len(f)):
                    high = dis[j + 1]
                    ai += MyTrap_tol(product(f[j], f_for_a), low, high, n=20, d=d)

                    low = high
                an.append(ai)

        elif method == 'simp':

            for i in range(1, N + 1):
                # i th coefficient of fourier series
                ai = 0

                f_for_a = lambda x: np.cos(i * np.pi * x / L)

                low = dis[0]

                for j in range(len(f)):
                    high = dis[j + 1]
                    ai += MySimp_tol(product(f[j], f_for_a), low, high, n=20, d=d)

                    low = high
                an.append(ai)

        return an

    def coef_bn():
 

        bn = []
        if method == 'quad':

            for i in range(1, N + 1):
                # i th coefficient of fourier series
                bi = 0

                f_for_b = lambda x: np.sin(i * np.pi * x / L)

                low = dis[0]

                for j in range(len(f)):
                    high = dis[j + 1]
                    bi += MyLegQuadrature_tol_1(product(f[j], f_for_b), low, high, n=5, d=d)

                    low = high

                bn.append(bi)

        elif method == 'trap':

            for i in range(1, N + 1):
                # i th coefficient of fourier series
                bi = 0

                f_for_b = lambda x: np.sin(i * np.pi * x / L)

                low = dis[0]

                for j in range(len(f)):
                    high = dis[j + 1]
                    bi += MyTrap_tol(product(f[j], f_for_b), low, high, n=20, d=d)

                    low = high
                bn.append(bi)

        elif method == 'simp':

            for i in range(1, N + 1):
                # i th coefficient of fourier series
                bi = 0

                f_for_b = lambda x: np.sin(i * np.pi * x / L)

                low = dis[0]

                for j in range(len(f)):
                    high = dis[j + 1]
                    bi += MySimp_tol(product(f[j], f_for_b), low, high, n=20, d=d)

                    low = high
                bn.append(bi)

        return bn

    if var == 0:
        a0 = coef_a0()
        an = coef_an()
        bn = np.zeros(N)

    elif var == 1:
        a0 = 0
        an = np.zeros(N)
        bn = coef_bn()

    else:
        a0 = coef_a0()
        an = coef_an()
        bn = coef_bn()

    return a0, an, bn


print('\nEnter the number of terms for the given functions:')
n = int(input())


def Q3i(n):
    f1 = [lambda x: 0, lambda x: 1]
    dis1 = [-1, 0, 1]
    L = 1
    a01, an1, bn1 = FourierCoeff("quad", 5, f1, dis1, 1, n, -1)
    T = [t for t in range(1, n + 1)]

    df = pd.DataFrame({'N': T, 'an1': an1, 'bn1': bn1})
    print('a01    : ', a01)
    print(df)
    S = []
    i = [1, 2, 5, 10, 20]
    xi = np.linspace(-2.5, 2.5, 50)

    y1 = []
    for n in i:
        sum_i = lambda x: a01 / 2

        for j in range(len(i)):
            f_add = lambda x: an1[j] * np.cos((j + 1) * x * np.pi / L) + bn1[j] * np.sin((j+1) * x * np.pi / L)

            sum_i = sumfn(sum_i, f_add)
        y1.append(sum_i)
        S.append(np.array(sum_i(xi)))

    s1 = []
    for k in range(len(S)):
        s1.append(sum(S[k]))
    dt = pd.DataFrame({'i': i, 'Sum_i': s1})
    print(dt)

    an_c = []
    bn_s = []

    def func(x, i):
        a0, an, bn = FourierCoeff("quad", d=5, f=f1, dis=dis1, L=1, N=i, var=-1)
        for j in range(0, len(an)):
            an_c.append(an[j] * np.cos((j + 1) * np.pi * x / 1))
            bn_s.append(bn[j] * np.sin((j + 1) * np.pi * x / 1))
        four = (a0 / 2 + sum(an_c) + sum(bn_s))
        return four, an_c, bn_s

    #plt.plot(xi, func(xi, 1)[0], label="i=1")
    #plt.plot(xi, func(xi, 2)[0], label="i=2")
    #plt.plot(xi, func(xi, 5)[0], label="i=5")
    #plt.plot(xi, func(xi, 10)[0], label="i=10")
    plt.plot(xi, func(xi, 20)[0], label="i=20")

    f = open("file1.dat", 'wb')
    x1 = np.column_stack((an1, bn1))
    np.savetxt(f, x1)
    f.close()

    x_true1 = [-2.5, -2, -1.92, -1, -0.92, 0, 0.08, 1, 1.08, 2, 2.08, 2.5]
    y_true1 = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]

    plt.plot(x_true1, y_true1, label='true_y')
    plt.legend()
    plt.show()


#print(Q3i(n))


def Q3ii(n):
    # for function (2)
    f2 = [lambda x: 0, lambda x: 1, lambda x: 0]
    dis2 = [-1, -0.5, 0.5, 1]
    a02, an2, bn2 = FourierCoeff("quad", 6, f=f2, dis=dis2, L=1, N=n, var=0)
    L = 1
    T = [t for t in range(1, n + 1)]
    xi = np.linspace(-2.5, 2.5, 50)

    df = pd.DataFrame({'N': T, 'an2': an2, 'bn2': bn2})
    print('a02 : ', a02)
    print(df)
    S = []
    i = [1, 2, 5, 10, 20]
    

    y1 = []
    for n in i:
        sum_i = lambda x: a02 / 2

        for j in range(len(i)):
            f_add = lambda x: an2[j] * np.cos((j + 1) * x * np.pi / L) + bn2[j] * np.sin((j+1) * x * np.pi / L)

            sum_i = sumfn(sum_i, f_add)
        y1.append(sum_i)
        S.append(np.array(sum_i(xi)))

    s1 = []
    for k in range(len(S)):
        s1.append(sum(S[k]))
    dt = pd.DataFrame({'i': i, 'Sum_i': s1})
    print(dt)

    an_c = []
    bn_s = []

    def func(x, i):
        a0, an, bn = FourierCoeff("quad", d=5, f=f2, dis=dis2, L=1, N=i, var=-1)
        for j in range(0, len(an)):
            an_c.append(an[j] * np.cos((j + 1) * np.pi * xi / 1))
            bn_s.append(bn[j] * np.sin((j + 1) * np.pi * xi / 1))
        four = (a0 / 2 + sum(an_c) + sum(bn_s))
        return four, an_c, bn_s

    plt.plot(xi, func(xi, 1)[0], label="i=1")
    #plt.plot(xi, func(xi, 2)[0], label="i=2")
    #plt.plot(xi, func(xi, 5)[0], label="i=5")
    #plt.plot(xi, func(xi, 10)[0], label="i=10")
    #plt.plot(xi, func(xi, 100)[0], label="i=20")
    #print(FourierCoeff("quad", d=5, f=f2, dis=dis2, L=1, N=10, var=-1))

    x_true1 = [-2.499999, -2.0000001, -1.500001, -1.4999999, -1.000001, -0.999999, -0.5000001, -0.49999999, 0.4999999, 0.5000001, 0.999999, 1.000001, 1.499999, 1.5000001, 2.0000001, 2.4999999]
    y_true1 = [1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1]
    f = open("file2.dat", 'wb')
    x1 = np.column_stack((an2, bn2))
    np.savetxt(f, x1)
    f.close()

    plt.plot(x_true1, y_true1, label='true_y')
    plt.legend()
    plt.show()
    #return a02, an2, bn2


#print(Q3ii(n))


def Q3iii(n):
    # for function(3)
    f3 = [lambda x: -0.5, lambda x: 0.5]
    dis3 = [-1, 0, 1]
    a03, an3, bn3 = FourierCoeff("quad", 6, f=f3, dis=dis3, L=1, N=n, var=1)
    L = 1
    T = [t for t in range(1, n + 1)]

    df = pd.DataFrame({'N': T, 'an3': an3, 'bn3': bn3})
    print('a03 : ', a03)
    print(df)
    S = []
    i = [1, 2, 5, 10, 20]
    xi = np.linspace(-2.5, 2.5, 50)

    y1 = []
    for n in i:
        sum_i = lambda x: a03 / 2

        for j in range(len(i)):
            f_add = lambda x: an3[j] * np.cos((j + 1) * x * np.pi / L) + bn3[j] * np.sin((j + 1) * x * np.pi / L)

            sum_i = sumfn(sum_i, f_add)
        y1.append(sum_i)
        S.append(np.array(sum_i(xi)))

    s1 = []
    for k in range(len(S)):
        s1.append(sum(S[k]))
    dt = pd.DataFrame({'i': i, 'Sum_i': s1})
    print(dt)

    an_c = []
    bn_s = []

    def func(x, i):
        a0, an, bn = FourierCoeff("quad", d=5, f=f3, dis=dis3, L=1, N=i, var=-1)
        for j in range(0, len(an)):
            an_c.append(an[j] * np.cos((j + 1) * np.pi * xi / 1))
            bn_s.append(bn[j] * np.sin((j + 1) * np.pi * xi / 1))
        four = (a0 / 2 + sum(an_c) + sum(bn_s))
        return four, an_c, bn_s

    #plt.plot(xi, func(xi, 1)[0], label="i=1")
    #plt.plot(xi, func(xi, 2)[0], label="i=2")
    #plt.plot(xi, func(xi, 5)[0], label="i=5")
    #plt.plot(xi, func(xi, 10)[0], label="i=10")
    #plt.plot(xi, func(xi, 20)[0], label="i=20")

    x_true1 = [-2.5, -2, -1.92, -1, -0.92, 0, 0.08, 1, 1.08, 2, 2.08, 2.5]
    y_true1 = [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5]

    f = open("file3.dat", 'wb')
    x1 = np.column_stack((an3, bn3))
    np.savetxt(f, x1)
    f.close()

    plt.plot(x_true1, y_true1, label='true_y')
    plt.legend()
    plt.show()
    #return a03, an3, bn3


print(Q3iii(n))



