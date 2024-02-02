import numpy as np
import matplotlib.pyplot as plt
from IVP_Module import *

def func_x(x,y_vec):
    ans_vec=np.zeros((2))
    ans_vec[0]=y_vec[1]
    ans_vec[1]=2*(y_vec[0]**3)
    return ans_vec

def row_cutting(mat,row_no):
    col_no=int(np.size(mat)/len(mat))
    new_mat=np.zeros((row_no,col_no))
    for i in range(row_no):
        new_mat[i,:]=mat[i,:]
    return new_mat

def analy_y(x):
    return 1/(x+3)

def analy_deriv_y(x):
    return -1/((x+3)**2)


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
        ans_mat=RK_4(func_x,para,t)
        last_val=a3*ans_mat[-1,0]+a4*ans_mat[-1,1]
        return abs(b2-last_val),ans_mat,t
    
    if tol==None:
        tol=-999
    
    s_k=[]
    y_mat_s=np.zeros((N_max,no_pt))
    y_mat_d_s=np.zeros((N_max,no_pt))
    
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
                
            
a=0
b=1
#part1 Dirichlet
alpha1=1
alpha2=0
beta1=1/3
alpha3=1
alpha4=0
beta2=1/4
s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,8,10**-6)
graph_sketching(t,y_mat_s,s_k,"x","y(x)","Dirichlet Conditions:y(0)=0.333 y(1)=0.25",analy_y)
graph_sketching(t,y_d_s,s_k,"x","y'(x)","Dirichlet Conditions:y(0)=0.333 y(1)=0.25",analy_deriv_y)


#part-2 Neumann Conditions
alpha1=0
alpha2=1
beta1=-1/9
alpha3=0
alpha4=1
beta2=-1/16
s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,8,10**-6,s_0=0.6,s_1=0.8)
graph_sketching(t,y_mat_s,s_k,"x","y(x)","Neumann Conditions:y'(0)=0.1111,y'(1)=-1/16",analy_y)
graph_sketching(t,y_d_s,s_k,"x","y'(x)","Neumann Conditions:y'(0)=0.1111,y'(1)=-1/16",analy_deriv_y)


#part-3 mixed with robin
alpha1=3
alpha2=-9
beta1=2
alpha3=1
alpha4=0
beta2=1/4
s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,8,10**-6)
graph_sketching(t,y_mat_s,s_k,"x","y(x)","robin a part",analy_y)
graph_sketching(t,y_d_s,s_k,"x","y'(x)","robin a part",analy_deriv_y)

#part-4 mixed with robin
alpha1=1
alpha2=0
beta1=1/3
alpha3=2
alpha4=2
beta2=3/48
s_k,ans_mat,y_mat_s,t,y_d_s=linear_shooting(a,b,alpha1,alpha2,alpha3,alpha4,beta1,beta2,func_x,8,10**-6)
graph_sketching(t,y_mat_s,s_k,"x","y(x)","robin b part",analy_y)
graph_sketching(t,y_d_s,s_k,"x","y'(x)","robin b part",analy_deriv_y)


