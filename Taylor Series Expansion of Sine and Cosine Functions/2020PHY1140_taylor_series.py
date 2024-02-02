#Name Preetpal Singh
#Roll No. 2020PHY1140
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
'''In MySinFunction user defined function, sin function is approximated with taylor expansion'''
def MySinFunction(x,n):
    y = []
    sol = 0
    for i in range(0,n,1):
        sol = sol + ((-1)**i * (x)**(2*i + 1))/math.factorial(2*i + 1)
        #y.append(sol)
    return sol
'''In MyCosFunction user defined function, cos function is approximated with taylor expansion'''    
def MyCosFunction(x,n):
    y = []
    sol = 0
    for i in range(0,n,1):
        sol = sol + ((-1)**i * (x)**(2*i))/math.factorial(2*i)
        #y.append(sol)
    return sol
'''This user defined function take domain and tolerance as input and return the output within tolerance limit and number of terms of expansion'''
def my_fun(x,tol):
    fun_cal,n,lis = 0,0,[]
    while True:  #Loop iterates till we get output within tolerance
        fun_cal = fun_cal+(-1)**n*x**((2*n)+1)/math.factorial((2*n)+1)
        lis.append(fun_cal)
        n+=1
        if len(lis)>=2: #list should have atleast two elements to compare
            if lis[-2]==0: #if 2nd last element of list is zero then break
                break
            else:
                err = abs((lis[-1]-lis[-2])/lis[-2]) #Calculate relative error
            if err <= tol: #if relative error is within tolerance then break
                break
    return fun_cal,n

'''This is another function in which tolerance is taken as input from user and output from my_func is used to print the table that compares results with inbuilt function with number of terms of expansion'''
def pr_values(x_a):
    tol = float(input("Input tolerance value:"))
    fun_cal=[];n=[]
    for x in x_a:
        h=my_fun(x,tol)
        fun_cal.append(h[0]);n.append(h[1])
    return fun_cal,n

'''graph function is to automate plotting for sin and cosine function. This user defined function takes arguments xrange(domain of function), inbuilt function(sinx,cosx) to compare, number of terms to expand taylor series, function(either MyCosFunction/MySinFunction), title, ylabel and title of subplot 1 and subplot 2 respectively'''
def graph(x_0,fun,m,req_func,title,ylabel,sub1_title,sub2_title,filename):
    fig, (ax1, ax2) = plt.subplots(1,2) 
    fig.suptitle(title, fontsize=20)
    ax1.plot(x_0,sin_x_0,color='green',label="Inbuilt")    
    for j in m:
            apr_sin = [req_func(ang,j) for ang in x_0]
            ax1.plot(x_0,apr_sin,label="n ={}".format(j),marker='o',markersize=3)
    ax1.legend();ax1.set_ylim([-5,5]);ax1.set_xlabel("x");ax1.grid()
    ax1.set_ylabel("{}".format(ylabel));ax1.set_title("{}".format(sub1_title))  
    m=np.arange(2,21,2)
    ax2.plot(m,np.sin(const_value(m)),color='green',label = "Inbuilt")

    for j in m:
            apr_sin = MySinFunction(const_value(j),j)
            ax2.plot(j,apr_sin,label="n ={}".format(j),marker='*',markersize=6)
    ax2.legend();ax2.set_xlabel("number of terms of taylor expansion");ax2.set_ylabel("{}".format(ylabel))
    ax2.grid();ax2.set_title("{}".format(sub2_title));plt.savefig("{}".format(filename))
    plt.show()

'''In const_value function @np.vectorize is used to return the same value if we iterate the function in loop'''
@np.vectorize
def const_value(m):    
    return np.pi/4

if __name__ == "__main__":
    '''--------------------------------------Q1-----------------------------------------------------'''
    x_0=np.linspace(-2*np.pi,2*np.pi,20) #Domain to calculate 
    sin_x_0 = np.sin(x_0)
    m = [1,2,5,10,20]
    graph(x_0,np.sin(x_0),m,MySinFunction,"Taylor Series Expansion of sin(x)","sin(x)","sin(x)","sin(x) at π/4","/home/hinton/Semester_4/MP3/Practical/Programmes/sin.jpg")
    graph(x_0,np.cos(x_0),m,MyCosFunction,"Taylor Series Expansion of cos(x)","cos(x)","cos(x)","cos(x) at π/4","//home/hinton/Semester_4/MP3/Practical/Programmes/cos.jpg")

    '''---------------------------------------Q2-----------------------------------------------'''
    x = np.arange(np.pi/8,np.pi,np.pi/8) 
    sin_in = np.sin(x)
    output = pr_values(x)
    data = {"x":x,"sin(x) cal":output[0],"n":output[1], "sin_inbuilt":sin_in}
    print(pd.DataFrame(data)) #print pandas dataframe 
    plt.scatter(x,output[0],c = "r",label="Value obtained by Taylor Expansion")
    plt.plot(x,sin_in,label="Inbuilt Function");plt.savefig("/home/hinton/Semester_4/MP3/Practical/Programmes/final.jpg") #save the output in image file
    plt.xlabel("x");plt.ylabel("sin(x)");plt.title("Taylor Expansion of Sine Function");plt.grid(True);plt.legend()
    plt.show()