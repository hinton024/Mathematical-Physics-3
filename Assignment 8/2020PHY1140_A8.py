from IVP import euler, RK_2, RK_4
import numpy as np
import matplotlib.pyplot as plt
from  scipy.integrate import RK45
from prettytable import PrettyTable
#Funcion To Be Defined(Not to be included in Module)
def func1(x,x_vec):
    ans_vec = np.zeros((3))
    ans_vec[0] = x_vec[1] - x_vec[2] + x
    ans_vec[1] = 3*x**2
    ans_vec[2] = x_vec[1] + np.exp(-x)
    return ans_vec
def graph(x,analytic,euler_final,rk2_final,rk4_final,title):
    fig,axs=plt.subplots(3,2,figsize=(15,15))
    fig.suptitle(title, fontsize=15)
    ax11,ax12,ax21,ax22,ax31,ax32=axs[0][0],axs[0][1],axs[1][0],axs[1][1],axs[2][0],axs[2][1]
    ax11.plot(x[0],euler_final[0],'^', color='green',label="euler"),ax11.plot(x[0],rk2_final[0],'-', color='black',label="rk2")
    ax11.plot(x[0],rk4_final[0],'>',color='black',label="rk4"),ax11.plot(x[0],analytic[0],'*',color='brown',label="Analytic")
    ax11.set_title("$N= 10^1$"),ax11.set_ylabel("Y"),ax11.set_xlabel("x")      
    ax12.plot(x[1],euler_final[1],'^', color='green',label="euler"),ax12.plot(x[1],rk2_final[1],'-', color='black',label="rk2")
    ax12.plot(x[1],rk4_final[1],'>',color='black',label="rk4"),ax12.plot(x[1],analytic[1],'*',color='brown',label="Analytic")
    ax12.set_title("$N=10^2$"),ax12.set_ylabel("y"),ax12.set_xlabel("x")      
    ax21.plot(x[2],euler_final[2],'^', color='green',label="euler"),ax21.plot(x[2],rk2_final[2],'-', color='black',label="rk2")
    ax21.plot(x[2],rk4_final[2],'>',color='black',label="rk4"),ax21.plot(x[2],analytic[2],'*',color='brown',label="Analytic")
    ax21.set_title("$N=10^3$"),ax21.set_ylabel("y"),ax21.set_xlabel("x")      
    ax22.plot(x[3],euler_final[3],'^', color='green',label="euler"),ax22.plot(x[3],rk2_final[3],'-', color='black',label="rk2")
    ax22.plot(x[3],rk4_final[3],'>',color='black',label="rk4"),ax22.plot(x[3],analytic[3],'*',color='brown',label="Analytic")
    ax22.set_title("$N=10^4$"),ax22.set_ylabel("y"),ax22.set_xlabel("x")      
    ax31.plot(x[4],euler_final[4],'^', color='green',label="euler"),ax31.plot(x[4],rk2_final[4],'-', color='black',label="rk2")
    ax31.plot(x[4],rk4_final[4],'>',color='black',label="rk4"),ax31.plot(x[4],analytic[4],'*',color='brown',label="Analytic")
    ax31.set_title("$N=10^5$"),ax31.set_ylabel("y"),ax31.set_xlabel("x")      
    ax32.plot(x[5],euler_final[5],'^', color='green',label="euler"),ax32.plot(x[5],rk2_final[5],'-', color='black',label="rk2")
    ax32.plot(x[5],rk4_final[5],'>',color='black',label="rk4"),ax32.plot(x[5],analytic[5],'*',color='brown',label="Analytic")
    ax32.set_title("$N=10^6$"),ax32.set_ylabel("y"),ax32.set_xlabel("x")      
    ax11.legend(),ax11.grid(True),ax12.legend(),ax12.grid(True),ax21.legend(),ax21.grid(True),ax22.legend(),ax22.grid(True)
    ax31.legend(),ax31.grid(True),ax32.legend(),ax32.grid(True);plt.tight_layout()
    plt.show()
x=np.linspace(0,1,100)
initial_conds = [1,1,-1]
euler_final,rk2_final,rk4_final,analytic_final,x_final=[],[],[],[],[]
# for j in np.arange(1,11,1.5):
E1,E2,E3=[],[],[];p=[]
for i in np.arange(1,3,1):
    x = np.linspace(0,2.5,10**i)
    analytic = [-0.05*x**5+0.25*x**4+x+2-np.exp(-x),x**3+1,0.25*x**4+x-np.exp(-x)]
    euler_1 = euler(func1,initial_conds,x).T[1]
    rk2 = RK_2(func1,initial_conds,x).T[1]
    rk4= RK_4(func1,initial_conds,x).T[1]
    euler_final.append(euler_1)
    rk2_final.append(rk2)
    rk4_final.append(rk4)
    analytic_final.append(analytic[1])
    x_final.append(x)
    #print(RK_4(func1,initial_conds, x))
    euler_error = np.max(abs(analytic[1]-(euler(func1,initial_conds, x)).T[1]))
    rk2error = np.max(abs(analytic[1]-(RK_2(func1,initial_conds, x)).T[1]))
    rk4error = np.max(abs(analytic[1]-(RK_4(func1,initial_conds, x)).T[1]))
    E1.append(rk2error);E2.append(rk4error);E3.append(euler_error)
    
    p.append(10**i)
# plt.scatter(np.log10(p),np.log10(E3),label="Euler upto N ={}".format(10**i))
# plt.scatter(np.log10(p),np.log10(E1),label="RK2 upto N ={}".format(10**i))
# plt.scatter(np.log10(p),np.log10(E2),label="RK4 upto N ={}".format(10**i))
# plt.legend();plt.grid(True);plt.xlabel("$log10(N)$");plt.ylabel("$log10(Error)$");
# plt.title("Error Plot for Differential equation 2 for limit [0,2.5]");plt.tight_layout()
# plt.show()
# graph(x_final,analytic_final,euler_final,rk2_final,rk4_final,"Y3 Differential Equation")

def func2(x,x_vec):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = 2*x_vec[1]-2*x_vec[0]+np.exp(2*x)*np.sin(x)
    return ans_vec    
initial_conds = [-0.4,-0.6]
x=np.linspace(0,1,6)
print(x)
print(RK_2(func2,initial_conds,x))

