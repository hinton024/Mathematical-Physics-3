import numpy as np
import matplotlib.pyplot as plt
from IVP import *
from math import *
import pandas as pd
from scipy.stats import linregress



def func_x(x,y_vec):
    ans_vec=np.zeros((2))
    ans_vec[0]=y_vec[1]
    ans_vec[1]=np.sin(3*x)+(-1*y_vec[0])
    return ans_vec

def row_cutting(mat,row_no):
    col_no=int(np.size(mat)/len(mat))
    new_mat=np.zeros((row_no,col_no))
    for i in range(row_no):
        new_mat[i,:]=mat[i,:]
    return new_mat

def graph_log_sketch(X,Y):
    fig,ax=plt.subplots()
    plt.plot(X,Y)
    ax.set_xscale("log",base=2)
    ax.set_yscale("log")
    ax.set_title("log(E) vs log(n)")
    ax.set_xlabel("log(N)")
    ax.set_ylabel("log(E)")
    plt.grid()
    plt.show()

def analy_y(x):
    return (3/8)*np.sin(x)-np.cos(x)-(1/8)*np.sin(3*x)

def analy_deriv_y(x):
    return (3/8)*np.cos(x)+np.sin(x)-(3/8)*np.cos(3*x)


def graph_sketching(x,y_mat,s_k,x_axis,y_axis,title,analy_y):
    fig,ax=plt.subplots()
    n=len(y_mat)
    for i in range(n):
        plt.plot(x,y_mat[i,:],"--",label="s="+str('%.4f'%(s_k[i])))
        plt.scatter(x,y_mat[i,:])
    y_true=[]
    analy_y=np.vectorize(analy_y)
    y_true=analy_y(x)
    plt.plot(x,y_true,label="analytic values")
    plt.grid()
    plt.legend()
    ax.set_title(title)
    ax.set_xlabel(x_axis)
    ax.set_ylabel(y_axis)
    plt.plot()
    plt.show()
    


def linear_shooting(a,b,a1,a2,a3,a4,b1,b2,func_x,no_pt,tol=None,s_0=0,s_1=1,N_max=50):
    
    def phi(s):
        if a2==0:
            para=[b1/a1,s]
        else:
            deriv_s=(b1-(a1*s))/a2
            para=[s,deriv_s]
        t=np.linspace(a,b,no_pt,float)
        ans_mat=RK_fourth_vec(t,para,func_x)
        last_val=a3*ans_mat[-1,0]+a4*ans_mat[-1,1]
        return abs(b2-last_val),ans_mat,t
    
    if tol==None:
        tol=-999
        
    s_k=[]
    y_mat_s=np.zeros((53,no_pt))
    y_mat_d_s=np.zeros((53,no_pt))
    
    err,ans_mat,t=phi(s_0)
    s_k.append(s_0)
    y_mat_s[0,:]=ans_mat[:,0]
    y_mat_d_s[0,:]=ans_mat[:,1]
    
    if err<tol or N_max==1:
        y_mat_s=row_cutting(y_mat_s,len(s_k))
        y_mat_d_s=row_cutting(y_mat_d_s,len(s_k))
        return s_k,ans_mat,y_mat_s,t,y_mat_d_s
    
    else:
        err,ans_mat,t=phi(s_1)
        s_k.append(s_1)
        y_mat_s[1,:]=ans_mat[:,0]
        y_mat_d_s[1,:]=ans_mat[:,1]
        
        if err<tol or N_max==2:
            y_mat_s=row_cutting(y_mat_s,len(s_k))
            y_mat_d_s=row_cutting(y_mat_d_s,len(s_k))
            return s_k,ans_mat,y_mat_s,t,y_mat_d_s
            
        else:
            step=2
            while step<N_max:
                s_2 = s_0 - (s_1-s_0)*phi(s_0)[0]/( phi(s_1)[0] - phi(s_0)[0] ) 
                s_k.append(s_2)
                s_0 = s_1
                s_1 = s_2
                step = step + 1 
                diff,ans_mat,t=phi(s_2)
                y_mat_s[step-1,:]=ans_mat[:,0]
                y_mat_d_s[step-1,:]=ans_mat[:,1]
                if diff<tol:
                    y_mat_s=row_cutting(y_mat_s,len(s_k))
                    y_mat_d_s=row_cutting(y_mat_d_s,len(s_k))
                    return s_k,ans_mat,y_mat_s,t,y_mat_d_s
    if tol!=-999:  
        print("tolerance not reached")            
    y_mat_s=row_cutting(y_mat_s,len(s_k))
    y_mat_d_s=row_cutting(y_mat_d_s,len(s_k))
    return s_k,ans_mat,y_mat_s,t,y_mat_d_s

def calc_rms_error(y_vec1,y_vec2):
    sum_ele=0
    for i in range(len(y_vec1)):
        ele=(y_vec1[i]-y_vec2[i])**2
        sum_ele=sum_ele+ele
    ans=sqrt(sum_ele/len(y_vec1))
    return ans
                
                
a=0
b=(np.pi)/2
#part1 Dirichlet
alpha1=1
alpha2=1
beta1=-1
alpha3=0
alpha4=1
beta2=1
s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,8,tol=10**-6,s_0=0,s_1=1)
graph_sketching(t,y_mat_s,s_k,"x","y(x)","Dirichlet Conditions:y(0)=0.333 y(1)=0.25",analy_y)
graph_sketching(t,y_d_s,s_k,"x","y'(x)","Dirichlet Conditions:y(0)=0.333 y(1)=0.25",analy_deriv_y)

s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,5,10**-6)
analy_y=np.vectorize(analy_y)
corr_list=analy_y(t)
err_list=np.abs(ans_mat[:,0]-corr_list)
data_mat1=np.column_stack((t,ans_mat[:,0],corr_list,err_list))
print("Number of points=4")
pdf1=pd.DataFrame(data_mat1,columns=["x_i","y_num","y_anal","error"])
print(pdf1)

print()

s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,9,10**-6)
analy_y=np.vectorize(analy_y)
corr_list=analy_y(t)
err_list=np.abs(ans_mat[:,0]-corr_list)
data_mat1=np.column_stack((t,ans_mat[:,0],corr_list,err_list))
print("Number of points=8")
pdf1=pd.DataFrame(data_mat1,columns=["x_i","y_num","y_anal","error"])
print(pdf1)



N_list=np.logspace(1.0,6.0,num=6,base=2)
err_abs=np.zeros(len(N_list))
err_rms=np.zeros(len(N_list))
ratio_err1=np.zeros(len(N_list))
ratio_err2=np.zeros(len(N_list))
marker_list=["s","o","d","+",".","v"]

count=0
for i in N_list:
    s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,int(i+1),tol=None,N_max=4)
    analy_y=np.vectorize(analy_y)
    err_abs[count]=np.max(np.abs(ans_mat[:,0]-analy_y(t)))
    err_rms[count]=calc_rms_error(ans_mat[:,0],analy_y(t))
    plt.scatter(t,ans_mat[:,0],label=str(i),marker=marker_list[count])
    if count==0:
        ratio_err1[count]=None
        ratio_err2[count]=None
    else:
        ratio_err1[count]=err_abs[count-1]/err_abs[count]
        ratio_err2[count]=err_rms[count-1]/err_rms[count]
    count=count+1
plt.plot(t,analy_y(t),label="exact solution")
data_mat2=np.column_stack((N_list,err_abs,ratio_err1,err_rms,ratio_err2))
print()
pdf2=pd.DataFrame(data_mat2,columns=["N","Absolute error","error ratio","RMS error","error ratio"])
print(pdf2)
plt.grid()
plt.legend()
plt.show()

graph_log_sketch(N_list,err_abs)
log_x=np.log10(N_list)
log_y=np.log10(err_abs)
print("slope,intercept:",linregress(log_x,log_y)[0:2])
