import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math as m
from scipy.stats import linregress


def finite_diff_method(a,b,alpha1,alpha2,alpha3,beta1,beta2,beta3,N,func_p,func_q,func_r):
    
    def func_ldu(x_arr):
        arr_l=np.array([])
        arr_d=np.array([])
        arr_u=np.array([])
        b_vec=np.array([])
        
        arr_l=-1-(h/2)*func_p(x_arr)
        arr_d=2+(h**2)*func_q(x_arr)
        arr_u=-1+(h/2)*func_p(x_arr)
        arr_r=func_r(x_arr)
        b_vec=(-1*(h**2)*arr_r)
        
        if alpha2==0:
            a11=1
            a12=0
            a_n1=1
            a_n2=0
            b1=(alpha3/alpha1)
            b_n1=(beta3/beta1)
            
        else:
            a11=arr_d[0]+((2*h*arr_l[0]*alpha1)/alpha2)
            a12=-2
            a_n1=arr_d[-1]-((2*h*arr_u[-1]*beta1)/beta2)
            a_n2=-2
            b1=(-1*(h**2))*arr_r[0]+((2*h*arr_l[0]*alpha3)/alpha2)
            b_n1=(-1*(h**2))*arr_r[-1]-((2*h*arr_u[-1]*beta3)/beta2)
        
        arr_d[0]=a11
        arr_d[-1]=a_n1
        arr_u[0]=a12
        arr_l[-1]=a_n2
        b_vec[0]=b1
        b_vec[-1]=b_n1
        
        return arr_l[1:],arr_d,arr_u[:-1],b_vec
    
    
    def form_tri_matrix(li,di,ui,x_arr):
        N=len(di)
        
        A_mat=np.zeros((N,N))
        
        for i in range(N):
            A_mat[i][i]=di[i]
            if i<N-1:
                A_mat[i][i+1]=ui[i]
                A_mat[i+1][i]=li[i]
        return A_mat
    

    def thomas_algo(A_mat,b_vec):
        N=len(b_vec)
        a=np.zeros((N))
        b=np.zeros((N))
        c=np.zeros((N))
        d=np.zeros((N))
        
        d=b_vec
        for i in range(N):
            b[i]=A_mat[i][i]
            if i<(N-1) :
                a[i+1]=A_mat[i+1][i]
                c[i]=A_mat[i][i+1]
        
        cp = np.zeros(N) 
        dp = np.zeros(N) 
        X = np.zeros(N) 
        
        cp[0] = c[0]/b[0]  
        dp[0] = d[0]/b[0]
        
        for i in np.arange(1,(N),1):
            dnum = b[i] - a[i]*cp[i-1]
            cp[i] = c[i]/dnum
            dp[i] = (d[i]-a[i]*dp[i-1])/dnum
    
        # Perform Back Substitution
        X[(N-1)] = dp[N-1]  # Obtain last xn 

        for i in np.arange((N-2),-1,-1):  # use x[i+1] to obtain x[i]
            X[i] = (dp[i]) - (cp[i])*(X[i+1])
    
        return(X)
          
    
    h=(b-a)/(N+1)
    x_arr=np.linspace(a,b,N+2,float)
    li,di,ui,bv=func_ldu(x_arr)
    A_mat=form_tri_matrix(li,di,ui,x_arr)
    omega=thomas_algo(A_mat,bv)
    
    return x_arr,omega

def calc_rms_error(y_vec1,y_vec2):
    sum_ele=0
    for i in range(len(y_vec1)):
        ele=(y_vec1[i]-y_vec2[i])**2
        sum_ele=sum_ele+ele
    ans=m.sqrt(sum_ele/len(y_vec1))
    return ans    


#Q1
def func_p1(x):
    ans_arr=np.zeros(len(x))
    for i in range(len(x)):
        ans_arr[i]=0
    return ans_arr
    
def func_q1(x):
    ans_arr=np.zeros(len(x))
    for i in range(len(x)):
        ans_arr[i]=(np.pi**2)
    return ans_arr
    
def func_r1(x):
    ans_arr=np.zeros(len(x))
    for i in range(len(x)):
        ans_arr[i]=-2*(np.pi**2)*np.sin(np.pi*x[i])
    return ans_arr
    
def analytic_x1(x):
    return np.sin(np.pi*x)

x_vals11, approximation11 = finite_diff_method(0, 1, 1, 0, 0, 1, 0, 0, 3 ,func_p1, func_q1,func_r1)
y_anal11 = analytic_x1(x_vals11)
abs_error11 = np.abs(y_anal11 - approximation11)
err_mat11=np.column_stack((x_vals11,approximation11,y_anal11,abs_error11))
print("N=3")
data_11=pd.DataFrame(err_mat11,columns=["x_i","y_num","y_analytic","error"])
print(data_11)


x_vals12, approximation12 = finite_diff_method(0, 1, 1, 0, 0, 1, 0, 0, 8 ,func_p1, func_q1,func_r1)
y_anal12 = analytic_x1(x_vals12)
abs_error12 = np.abs(y_anal12 - approximation12)
err_mat12=np.column_stack((x_vals12,approximation12,y_anal12,abs_error12))
print("N=8")
data_12=pd.DataFrame(err_mat12,columns=["x_i","y_num","y_analytic","error"])
print(data_12) 

N1 = []
max_abs_error = []
max_abs_error_ratio = [0]
rms_error = []
rms_error_ratio = [0]

for i in range(1, 7):
    N = 2**i
    x_vals1, approximation1 = finite_diff_method(0, 1, 1, 0, 0, 1, 0, 0, N ,func_p1, func_q1,func_r1)
    y_anal1 = analytic_x1(x_vals1)
    abs_error1 = np.abs(y_anal1 - approximation1)
    rms_error1 = calc_rms_error(y_anal1, approximation1)
    max_abs_error1 = np.max(abs_error1)
    plt.plot(x_vals1, approximation1, label = "N={}".format(N), linestyle='dashed')
    plt.scatter(x_vals1, approximation1, s = 10)

    N1.append(N)
    max_abs_error.append(max_abs_error1)
    rms_error.append(rms_error1)

x = np.linspace(0, 1, 100)
plt.plot(x, analytic_x1(x), label = 'Analytic Solution')
plt.title('Variation of solution with N')
plt.xlabel('x')
plt.ylabel('Solution(y)')
plt.legend()
plt.grid()
plt.show()

for i in range(0,5):
    ratio1 = max_abs_error[i]/max_abs_error[i+1]
    max_abs_error_ratio.append(ratio1)

    ratio2 = rms_error[i]/rms_error[i+1]
    rms_error_ratio.append(ratio2)


convergence_data1 =np.column_stack((N1,max_abs_error ,max_abs_error_ratio ,rms_error,rms_error_ratio))
convergence_table1=pd.DataFrame(convergence_data1,columns=["N","max_abs_error","Error Ratio","Rms Error","Error Ratio"])
print(convergence_table1)

plt.plot(N1, max_abs_error, label = 'Error')
plt.scatter(N1, max_abs_error)
plt.xscale('log')
plt.yscale('log')
plt.title('Log Plot')
plt.xlabel('N')
plt.ylabel('Max absolute error')
plt.legend()
plt.grid()
plt.show()

log_x=np.log10(N1)
log_y=np.log10(max_abs_error)
print("slope,intercept:",linregress(log_x,log_y)[0:2])


#Q2

def func_p2(x):
    ans_arr=np.zeros(len(x))
    for i in range(len(x)):
        ans_arr[i]=0
    return ans_arr
    
def func_q2(x):
    ans_arr=np.zeros(len(x))
    for i in range(len(x)):
        ans_arr[i]=-1
    return ans_arr
    
def func_r2(x):
    ans_arr=np.zeros(len(x))
    for i in range(len(x)):
        ans_arr[i]=np.sin(3*x[i])
    return ans_arr
    
def analytic_x2(x):
    return (3/8)*np.sin(x)-np.cos(x)-(1/8)*np.sin(3*x)


x_vals21, approximation21 = finite_diff_method(0, np.pi/2, 1, 1, -1, 0, 1, 1, 3 ,func_p2, func_q2,func_r2)
y_anal21 = analytic_x2(x_vals21)
abs_error21 = np.abs(y_anal21 - approximation21)
err_mat21=np.column_stack((x_vals21,approximation21,y_anal21,abs_error21))
print("N=3")
data_21=pd.DataFrame(err_mat21,columns=["x_i","y_num","y_analytic","error"])
print(data_21)


x_vals22, approximation22 = finite_diff_method(0, np.pi/2, 1, 1, -1, 0, 1, 1, 8 ,func_p2, func_q2,func_r2)
y_anal22 = analytic_x2(x_vals22)
abs_error22 = np.abs(y_anal22 - approximation22)
err_mat22=np.column_stack((x_vals22,approximation22,y_anal22,abs_error22))
print("N=8")
data_22=pd.DataFrame(err_mat22,columns=["x_i","y_num","y_analytic","error"])
print(data_22)

N2 = []
max_abs_error = []
max_abs_error_ratio = [0]
rms_error = []
rms_error_ratio = [0]

for i in range(1, 7):
    N = 2**i
    x_vals2, approximation2 = finite_diff_method(0, np.pi/2, 1, 1, -1, 0, 1, 1, N ,func_p2, func_q2,func_r2)
    y_anal2 = analytic_x2(x_vals2)
    abs_error2 = np.abs(y_anal2 - approximation2)
    rms_error2 = calc_rms_error(y_anal2, approximation2)
    max_abs_error2 = np.max(abs_error2)
    plt.plot(x_vals2, approximation2, label = "N={}".format(N), linestyle='dashed')
    plt.scatter(x_vals2, approximation2, s = 10)

    N2.append(N)
    max_abs_error.append(max_abs_error2)
    rms_error.append(rms_error2)

x = np.linspace(0, np.pi/2, 100)
plt.plot(x, analytic_x2(x), label = 'Analytic Solution')
plt.title('Variation of solution with N')
plt.xlabel('x')
plt.ylabel('Solution(y)')
plt.legend()
plt.grid()
plt.show()


for i in range(0,5):
    ratio1 = max_abs_error[i]/max_abs_error[i+1]
    max_abs_error_ratio.append(ratio1)

    ratio2 = rms_error[i]/rms_error[i+1]
    rms_error_ratio.append(ratio2)


convergence_data2 =np.column_stack((N2,max_abs_error ,max_abs_error_ratio ,rms_error,rms_error_ratio))
convergence_table2=pd.DataFrame(convergence_data2,columns=["N","max_abs_error","Error Ratio","Rms Error","Error Ratio"])
print(convergence_table2)

plt.plot(N2, max_abs_error, label = 'Error')
plt.scatter(N2, max_abs_error)
plt.xscale('log')
plt.yscale('log')
plt.title('Log Plot')
plt.xlabel('N')
plt.ylabel('Max absolute error')
plt.legend()
plt.grid()
plt.show()

log_x=np.log10(N2)
log_y=np.log10(max_abs_error)
print("slope,intercept:",linregress(log_x,log_y)[0:2])