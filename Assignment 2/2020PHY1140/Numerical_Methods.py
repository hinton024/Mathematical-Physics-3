from Numerical_Methods import My_Trap,My_Simp,MyLegQuadrature, MyLegQuadrature_tol,MySimp_tol,MyTrap_tol,MyLegQuadrature_tol_1
import numpy as np
import math
import matplotlib.pyplot as plt
import texttable as tt
tab = tt.Texttable()
# 3 b i
#MyTrap_tol(My_Trap, 1, 10, 100, 8)
#print("Here, Trapezoidal method in limit [1,10] gives exact results for linear function(x)")


#MyTrap_tol(My_Trap, 1, 10, 100, 8)
#print("Here, Trapezoidal method in [1,10] limit doesn't give exact result for polynomial of order 2")

# 3 b ii

#MySimp_tol(My_Simp, 1, 10, 100, 8)
#print("Here, Simpson method in limit [1,10] gives exact results for polynomial(x) of degree 3")


#MySimp_tol(My_Simp, 1, 10, 100, 8)
#print("Here, Simpson method in [1,10] limit doesn't give exact result for polynomial of order 4")


# 3 b iii

#MyLegQuadrature_tol("x", 1, 2, 2, 1, 8)
#print("Here, for m = 1 composite 2 point method of gauss quadrature in limit [1, 2] is accurate for polynomial of degree 3")


#MyLegQuadrature_tol("x", 1, 2, 2, 1, 8)
#print("Here, for m = 1 composite 2 point method of gauss quadrature in limit [1, 2] is not accurate for polynomial of degree 4")


#MyLegQuadrature_tol("x", 1, 2, 4, 1, 8)
#print("Here, for m = 1 composite 4 point method of gauss quadrature in limit [1, 2] is accurate for polynomial of degree 7")

#MyLegQuadrature_tol("x", 1, 2, 4, 1, 8)
#print("Here, for m = 1 composite 4 point method of gauss quadrature in limit [1, 2] is not accurate for polynomial of degree 8")


# 3 c
#Trapezoidal method table
m = np.arange(1,17,1)
n = 2*m
h = (1-0)/n
my_pi_1 = []
my_pi_2 = []
const = []
err_1, err_2 = [], []
for i in n:
    p = My_Simp("x", 0, 1, i)
    q = My_Trap("x", 0, 1, i)
    const.append(np.pi)
    e_1 = abs(p - np.pi)/np.pi
    e_2 = abs(q-np.pi)/np.pi
    err_1.append(e_1)
    err_2.append(e_2)
    my_pi_1.append(p)
    my_pi_2.append(q)
my_pi = [my_pi_1, my_pi_2]
#print(my_pi)
#plt.plot(n, err_1)
#plt.plot(np.log(h), np.log(err_1), c='r')
#plt.plot(np.log(h), np.log(err_2), c='b')
#plt.plot(n, err_2, c='b')
#plt.plot(n, const, c='b')
plt.grid(True)
plt.xlabel("log(h)")
plt.ylabel('log(error)')
plt.title("log plot for error vs step size")
plt.legend(["simpson", "trapezoidal"])
#plt.show()

# 3 d
'''headings_1 = ["my_pin(n)(simpson)","n","E = |my_pi(n)-pi|/pi"]
tab.header(headings_1)
for row in zip(my_pi_1, n, err_1):
    tab.add_row(row)
    tab.set_max_width(0)
    tab.set_precision(6)
s = tab.draw()
print(s)
tab.reset()'''
#mysimpson method table
'''headings_2 = ["my_pin(n)(trapezoidal)","n","E = |my_pi(n)-pi|/pi"]
tab.header(headings_2)
for row in zip(my_pi_2,n,err_2):
    tab.add_row(row)
    tab.set_max_width(0)
    tab.set_precision(6)
s = tab.draw()
print(s)
tab.reset()'''

# 3 e
n_ = [2, 4, 8, 16, 32, 64]
m_ = [1, 2, 4, 8, 16, 32]
matrix_pi_quad = []
e = []
'''for i in n_:
    for j in m_:
        u = MyLegQuadrature("x", 0, 1, i, j)
        matrix_pi_quad.append(u)
        error = abs(u - np.pi)/np.pi
        e.append(error)'''
for i in m_:
    for j in n_:
        u = MyLegQuadrature("x", 0, 1, i, j)
        matrix_pi_quad.append(u)
        error = abs(u - np.pi)/np.pi
        e.append(error)

t = np.reshape(matrix_pi_quad,(len(m_), len(n_)))
t_ = np.reshape(e,(len(m_), len(n_)))
#print(t)
#print(t_)

headings_3 = ["n","m=1","m=2","m=4","m=8","m=16","m=32"]
tab.header(headings_3)
for row in zip(n,t[0],t[1],t[2],t[3],t[4],t[5]):
    tab.add_row(row)
    tab.set_max_width(0)
    tab.set_precision(6)
s = tab.draw()
#print(s)
tab.reset()

# 3 e ii
'''plt.plot(n_, t_[0])
plt.plot(n_, t_[3])
plt.style.use('seaborn')
plt.grid(True)
plt.xlabel("n")
plt.ylabel('e(n)')
plt.title("e vs n")
plt.legend(["m = 1", "m = 8"])
plt.show()'''

plt.plot(m_, t_[0])
plt.plot(m_, t_[2])
plt.style.use('seaborn')
plt.grid(True)
plt.xlabel("m")
plt.ylabel('e(m)')
plt.title("e vs m")
plt.legend(["n = 2", "n = 8"])
plt.show()


'''plt.style.use('seaborn')
fig1, ax1 = plt.subplots()
ax1.plot(n_, t[0], color='#444444',
         linestyle='--', label='m=1')
ax1.plot(n_,t[3], label='m=8')
ax1.set_title('n vs pi')
ax1.set_xlabel('n')
ax1.set_ylabel('pi')
ax1.legend()
plt.tight_layout()'''
#plt.show()


#@np.vectorize
'''def const_value(m):
    return np.pi
for i in m:
    # cons_pi = []
    cons_pi = const_value(m)
plt.style.use('seaborn')
fig1, ax1 = plt.subplots()
ax1.plot(m_, t[1], color='#444444',
         linestyle='--', label='n=2')
ax1.plot(m_,t[3], label='n=8')
ax1.plot(m,cons_pi,label = "inbuilt method")
ax1.set_title('m vs pi')
ax1.set_xlabel('m')
ax1.set_ylabel('pi')
ax1.legend()
plt.tight_layout()
plt.show()'''


n_ = [2, 4, 8, 16, 32]
d = [1, 2, 3, 4, 5, 6, 7, 8]
matrix_pi_quad = []
for i in d:
    for j in n:
        u = MyLegQuadrature_tol_1("x", 0, 1, j, i)
        t = matrix_pi_quad.append(u)



g = np.reshape(t, (len(n_), len(d)))
print(g)
print(matrix_pi_quad)


