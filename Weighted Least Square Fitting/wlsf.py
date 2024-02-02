import matplotlib.pyplot as plt
import numpy as np

File_data = np.loadtxt("/home/hinton/Semester_4/MP3/Practical/p/1140.dat",skiprows=1,usecols=[1,2,3,4,5,6,7,8,9,10,11],delimiter=',',dtype='float')
p = [];q = [];r=[];s_dev=[];
for i in np.arange(0,10):
    mean = ((np.average((File_data[i][1:np.shape(File_data)[1]]))))
    standard_deviation = ((np.sum((File_data[i][1:np.shape(File_data)[1]]-mean)**2)/90))*4*(mean)**2
    weight =(1)/(standard_deviation)
    q.append(mean**2)
    r.append(File_data[i][0])
    p.append(weight)
    s_dev.append(standard_deviation)
print("masses",r)
print("y_mean",q)
print("weight",p)
N = np.array([r,q,p])
np.savetxt('/home/hinton/Semester_4/MP3/Practical/p/1140.csv',N,delimiter=',')

def Mywlsf(x_i,y_i,w_i):    
    Dell = (((np.sum(w_i))*(np.sum(w_i*x_i*x_i)))-(np.sum(w_i*x_i))**2)
    c = (((np.sum(w_i*(x_i**2)))*(np.sum(w_i*y_i)))-((np.sum(w_i*x_i))*(np.sum(w_i*y_i*x_i))))/Dell
    m = (((np.sum(w_i))*(np.sum(w_i*y_i*x_i)))-((np.sum(w_i*x_i))*(np.sum(w_i*y_i))))/Dell
    Error_in_intercept = np.sqrt((np.sum(w_i*x_i*x_i))/Dell)
    Error_in_slope = np.sqrt((np.sum(w_i))/Dell)   
    chi_sq = np.sum(((y_i-(m*x_i+c))**2)*w_i)
    Y_cal = m*x_i+c
    print("intercept = ",c)
    print("slope = ",m)
    print("Error_in_intercept = ",Error_in_intercept)
    print("Error_in_slope = ",Error_in_slope)
    print("chi_sq = ",chi_sq)
    print("Errors = ",Y_cal-y_i)
    print("Sum of residuals = ",(np.sum((Y_cal-y_i)*w_i)))
    print("Sum of square of residuals", (np.sum(w_i*(Y_cal-y_i))**2))
    print("Fitted Values = ",Y_cal)
    print("---------------------------------------------------------")
    print("spring constant = ", 4*np.pi**2/m)
    print("m = ", (4*(np.pi**2)/m)*(c/(4*np.pi**2)))
    fig1, (ax1,ax2) = plt.subplots(2,1)
    ax1.scatter(x_i,y_i,label="observed points")
    ax1.plot(x_i,Y_cal,label="fitted points")
    ax1.set_ylabel('$T^2$')
    ax1.set_title("Wlsf")
    ax1.errorbar(x_i, y_i, (Y_cal-y_i),fmt='o',label="error")
    ax1.legend(loc='upper left')
    ax1.grid(True)
    ax2.scatter(x_i,y_i,label="observed points")
    ax2.plot(x_i,Y_cal,label="fitted points")
    ax2.set_ylabel('$T^2$')
    ax2.set_xlabel('$M$')
    ax2.set_ylim(0.57,0.59)
    ax2.set_xlim(154,156)
    ax2.errorbar(x_i, y_i, (Y_cal-y_i),fmt='o',label="error")
    ax2.legend(loc='upper left')
    ax2.grid(True)
    plt.savefig("/home/hinton/Semester_4/MP3/Practical/p/1140.png")
    plt.show()
    
def Mylsf(x_i,y_i):
    print("--------------------------------------------------------")   
    print("Least Square Fitting") 
    Dell = np.size(x_i)*(np.sum(x_i**2))-(np.sum(x_i))**2
    c = (((np.sum(x_i**2))*(np.sum(y_i)))-((np.sum(x_i))*(np.sum(y_i*x_i))))/Dell
    m = (np.size(x_i)*np.sum(x_i*y_i)-np.sum(x_i)*np.sum(y_i))/Dell
    Y_cal = m*x_i+c
    S = np.sqrt(np.sum((y_i-Y_cal)**2)/(np.size(x_i)-2))
    Error_in_intercept = S*np.sqrt((np.sum(x_i**2)/(np.size(x_i)*np.sum(x_i**2)-(np.sum(x_i))**2)))
    Error_in_slope = S*np.sqrt((np.size(x_i)/(np.size(x_i)*np.sum(x_i**2)-(np.sum(x_i)**2))))
    chi_sq = np.sum(((y_i-(m*x_i+c))**2))
    Y_cal = m*x_i+c
    print("intercept = ",c)
    print("slope = ",m)
    print("Error_in_intercept = ",Error_in_intercept)
    print("Error_in_slope =",Error_in_slope)
    print("chi_sq = ",chi_sq)
    print("Errors = ",Y_cal-y_i)
    print("Sum of residuals = ",np.sum(y_i-(m*x_i+c)))
    print("Sum of square of residuals", np.sum(((y_i-(m*x_i+c))**2)))
    print("Fitted Values =",Y_cal)
    fig1, (ax1,ax2) = plt.subplots(2,1)
    ax1.scatter(x_i,y_i,label="observed points")
    ax1.plot(x_i,Y_cal,label="fitted points")
    ax1.set_ylabel('$T^2$')
    ax1.set_title("LSF")
    ax1.errorbar(x_i, y_i, Y_cal-y_i,fmt='o',label="error")
    ax1.legend(loc='upper left')
    ax1.grid(True)
    ax2.scatter(x_i,y_i,label="observed points")
    ax2.plot(x_i,Y_cal,label="fitted in y")
    ax2.set_ylabel('$T^2$')
    ax2.set_xlabel('$M$')
    ax2.set_ylim(0.57,0.59)
    ax2.set_xlim(154,156)
    ax2.errorbar(x_i, y_i, (Y_cal-y_i),fmt='o',label="error")
    ax2.legend(loc='upper left')
    plt.grid(True)
    plt.savefig("/home/hinton/Semester_4/MP3/Practical/p/1140_1.png")
    plt.show()
    

Mywlsf(N[0],N[1],N[2])
Mylsf(N[0],N[1])
