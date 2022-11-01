

import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import os
import multiprocessing
import subprocess
from matplotlib.backends.backend_pdf import PdfPages


color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']



# Adding the bias noise

Data = None
path = '/home/varghese/Desktop/20220619/hpf'
order_list = sorted(os.listdir(path))

Len = len(order_list)

for directory in order_list:
    file_path = os.path.join(path,directory,'fits')
    file_list = sorted(os.listdir(file_path))
    First_frame = file_list[0]
    File = os.path.join(file_path,First_frame)

    if Data is None:
        Data = fits.getdata(File)
    else:
        Data = Data + fits.getdata(File)

Avg_Data = Data/Len    







def objective(x,l,m,n):
    return l*x**2+m*x+n

def exponential(t,a,b):
    return a*(np.exp(-b*t)-1)

def ab_from_pq(p,q):
    c = p+(3/8)*q
    d=q

    b=d/c
    a=c/b
    return a,b


def v1(x):
    return x

def v2(x):
    return ((-3/8)*x+x**2/2)


def polynomial(t,p,q):
                                                                                                                                                                  
    V1 = v1(t)
    V2 = v2(t)

    return -p*V1+q*V2


def Constants(l,m,n):  # In terms of p                                                                                                                               

    def exp_p(t,p):
        #q = (l*p**2+m*p+n)                                                                                                                                          
        #c = (p+(3/8)*(l*p**2+m*p+n))                                                                                                                                
        #b= ((l*p**2+m*p+n) / (p+(3/8)*(l*p**2+m*p+n)))                                                                                                              
        #a= ((p+(3/8)*(l*p**2+m*p+n))**2/ (l*p**2+m*p+n))                                                                                                            

        return ((p+(3/8)*(l*p**2+m*p+n))**2/ (l*p**2+m*p+n))*(np.exp(-((l*p**2+m*p+n) / (p+(3/8)*(l*p**2+m*p+n)))*t)-1)
    return exp_p

def Constants1(l,m,n): # In terms of q                                                                                                                               

    def exp_p(t,p):
        # c = (((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q)                                                                                                         
        # d = q                                                                                                                                                      
        # b = d/c = (q/(((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))                                                                                               
        # a = c**2 / d = ((((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))**2/(q)                                                                                     
        return ((((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))**2/(q)*(np.exp(-(q/(((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))*t)-1)
    return exp_p



def curve_from_p(t,p,l,m,n):
    q = l*p**2+m*p+n
    c = p+(3/8)*q
    b=q/c
    a= c**2/q

    return a*(np.exp(-b*t)-1)




def Fitting_1(path,name):
    
    pvs = []
    qvs = []
    er_p = []
    er_q = []
    counts = []
    for channel in range(4):
        
        Filename = name+'{}.npy'.format(channel)
        File = os.path.join(path,Filename)
        dc = np.load(File)

        Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
        Sigma1 = np.median(np.sqrt(Sigma2.flatten()))+np.zeros(20)
        fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:512*(channel+1)]

        interval = 250
        
        for flux_bin in range(500,60000,interval):
            Avg_Data_Channel = Avg_Data[:,channel*512:channel*512+512]
            pix_mask = (dc[-1] > (flux_bin-interval//2)) & (dc[-1] < (flux_bin+interval//2)) & (fiber_mask == 1)
            avg_count = np.nanmedian(dc[10,pix_mask])
            delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/(dc[10,pix_mask]+np.median(Avg_Data_Channel[pix_mask])),axis=1)
            time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / (np.size(delta_change_cube,axis=0))

            if not np.isnan(avg_count) and delta_change_cube[-1]<0:

                Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

                try:
                    fitting_data_variables, error = curve_fit(polynomial,time1,delta_change_cube,np.array([1.0,0.5]),sigma = Sigma)

                    p,q = fitting_data_variables

                    dp = np.sqrt(error[0][0])
                    dq = np.sqrt(error[1][1])

                    pvs.append(p)
                    qvs.append(q)
                    counts.append(flux_bin)
                    er_p.append(dp)
                    er_q.append(dq)

                except Exception as e:
                    print(e)

                    print('{0},{1}'.format(avg_count,channel),'cannot be done')



    print("channel", channel, 'done pvs and qvs')


    return pvs,qvs,counts,er_p,er_q




def Fitting_2(channel,path,name):
    
    pvs,qvs,counts,er_p,er_q = Fitting_1(path,name)


    
    fitted_values,error = curve_fit(objective,pvs,qvs,np.array([0.5,0.5,0.0]),sigma=er_q)

    l,m,w = fitted_values
    er_l = np.sqrt(error[0][0])
    er_m = np.sqrt(error[1][1])
    er_e=np.sqrt(error[2][2])
    print(channel,'l:',l,'m:',m,'w:',w,'errors',er_l,er_m,er_e)                                                                                                                                      
    #a_n, b_n = ab_from_pq(pvs,objective(pvs,l,m,w))
    # Here, p is q actually                                                                                                                                          
    p_values = []
    p_err = []
    Count_new = []
    Avg_Data_Channel = Avg_Data[:,channel*512:channel*512+512]
    Filename = name+'{}.npy'.format(channel)
    File = os.path.join(path,Filename)
    dc = np.load(File)
    
    pdfplots = PdfPages("Decay_curves_ch_{0}.pdf".format(channel))
    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)
    interval = 250
    flag = 0
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    for flux_bin in range(500,60000,interval):
        
        
        pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc[10,pix_mask])
        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/(dc[10,pix_mask]+np.median(Avg_Data_Channel[pix_mask])),axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1)
        time1=time1/time1[-1]
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            P_ini=0.05
            flag+=1
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

            try:
                PQ,er_pq = curve_fit(polynomial,time1,delta_change_cube,np.array([0.05,0.05]),sigma=Sigma)
                P_new,Q_new = PQ

                ab, er_ab = curve_fit(exponential,time1,delta_change_cube,np.array([0.1,0.01]),sigma=Sigma)
                a,b = ab

                P_values, error = curve_fit(Constants(l,m,w),time1,delta_change_cube,np.array([P_ini]),sigma=Sigma)
                P = P_values[0]
                perr =  np.sqrt(error[0][0])
#                print('channel',channel,'p found',P,'initial p',P_ini,'error',perr)
                p_values.append(P)
                p_err.append(perr)
                Count_new.append(flux_bin)

                
                fig = plt.figure(figsize=(16,9))
                plt.plot(time1,delta_change_cube,'o-',label='data')
                        
                plt.plot(time1,exponential(time1,a,b),'--',label='exponential from a and b')   # Polynomial with a and b                                            
                plt.plot(time1,curve_from_p(time1,P,l,m,w),'--',label='exp finding p={}'.format(P))
                #plt.plot(time1,curve_from_p(time1,P_ini,l,m,w),'--',label='exp model from initial value of p={}'.format(P_ini))                                    
                plt.plot(time1,polynomial(time1,P_new,Q_new),'--',label='polynomial with p={0},q={0}'.format(P_new,Q_new)) # Polynomial with p and q              
                plt.plot(time1,polynomial(time1,P_new,objective(P_new,l,m,w)),'--',label='polynomial with p only'.format(P_new)) # Polynomial with p only           
                plt.xlabel('time')
                plt.legend()
                plt.ylabel('count')
                        
                plt.title('p Fitting Channel:{0}, Count{1}'.format(channel,flux_bin))
                pdfplots.savefig()
                        

                
                  


            except Exception as e:
                print(e)
    print(channel,'done pvs')
    pdfplots.close()
    #PVs[channel] = p_values
    #Counts[channel] = Count_new
    #er_P[channel] = p_err
    np.save('new_count_{}.npy'.format(channel),Count_new)
    np.save('new_pvs_{}.npy'.format(channel),p_values)
    np.save('new_p_err_{}.npy'.format(channel),p_err)
    



#Plotting

path = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G'

File = 'AvgDataCube_0G_'


#DS5: AvgDataCube_0G_


#Plotting

process1 = multiprocessing.Process(target=Fitting_2,args = [0,path,File])
process2 = multiprocessing.Process(target=Fitting_2,args = [1,path,File])
process3 = multiprocessing.Process(target=Fitting_2,args = [2,path,File])
process4 = multiprocessing.Process(target=Fitting_2,args = [3,path,File])


process1.start()
process2.start()
process3.start()
process4.start()

process1.join()
process2.join()
process3.join()
process4.join()





# Plotting


color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']
Shape = ['.','<','*','s']
plt.figure()
for channel in range(4):
    counts = np.load('new_count_{}.npy'.format(channel))
    pvs = np.load('new_pvs_{}.npy'.format(channel))
    err = np.load('new_p_err_{}.npy'.format(channel))
    print(np.median(err))
    plt.errorbar(counts,pvs,yerr=err,fmt=Shape[channel],color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
plt.ylim(-0.05,0.6)
plt.legend()
plt.xlabel('counts')
plt.ylabel('p')
plt.title('p vs counts from exponential,DS5-19/5')
plt.savefig('DS5_p_from_exp.png')

subprocess.run('rm -r *.npy', shell=True)
