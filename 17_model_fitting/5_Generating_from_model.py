

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time
import multiprocessing

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

t1 = time.time()

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
#    print(t)                                                                                                                                                                                              
    V1 = v1(t)
    V2 = v2(t)

    return -p*V1+q*V2


def Fitting(channel):
    pvs=np.load('pvs/pvs_{}.npy'.format(channel))
    qvs=np.load('qvs/qvs_{}.npy'.format(channel))
    er_p=np.load('er_p/er_p_{}.npy'.format(channel))
    er_q=np.load('er_q/er_q_{}.npy'.format(channel))
    avs=np.load('avs/avs_{}.npy'.format(channel))
    bvs=np.load('bvs/bvs_{}.npy'.format(channel))
    er_a=np.load('er_a/er_a_{}.npy'.format(channel))
    er_b=np.load('er_b/er_b_{}.npy'.format(channel))

    counts = np.load('counts/counts_{}.npy'.format(channel))
    
    a_n, b_n = ab_from_pq(pvs,qvs)

    fitted_values,error = curve_fit(objective,pvs,qvs,np.array([0.5,0.5,0.0]),sigma=er_q)
    
    l,m,e = fitted_values
    er_l = np.sqrt(error[0][0])
    er_m = np.sqrt(error[1][1])
    er_e=np.sqrt(error[2][2])
    print(channel,l,m,e,error)
    a_n, b_n = ab_from_pq(pvs,objective(pvs,l,m,e))
    

    for n in range(len(counts)):
        c = counts[n]
        delta_change_cube=np.load('Ch_{0}/Decay_curve_{1}_{2}.npy'.format(channel,channel,c))

        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / (np.size(delta_change_cube,axis=0))

        plt.figure(figsize=(16,9))
        plt.plot(time1,delta_change_cube,label='Data')
        plt.plot(time1,exponential(time1,avs[n],bvs[n]),label='Exponential_fit')
        plt.plot(time1,polynomial(time1,pvs[n],qvs[n]),label='Polynomian fit')
        plt.plot(time1,polynomial(time1,pvs[n],objective(pvs[n],l,m,e)),label='Polynomial from model')
        
        plt.plot(time1,exponential(time1,a_n[n],b_n[n]),label='Generated from model')
        plt.legend()
        plt.title('count={6},a={0},b={1},p={2},q={3},a_n={4},b_n={5}'.format(avs[n],bvs[n],pvs[n],qvs[n],a_n[n],b_n[n],c))
        plt.savefig('Ch_{0}/Model_comparison_{1}_{2}.png'.format(channel,channel,c))
        
    
    #plt.figure()
    #plt.errorbar(pvs,qvs,xerr=er_p,yerr=er_q,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
    #plt.plot(pvs,objective(pvs,l,m,n),'--',color=color_list[channel])
    #plt.savefig("Fitted_model_{}.png".format(channel))
    #print("Saved")

    
    

process1 = multiprocessing.Process(target=Fitting,args = [0])
process2 = multiprocessing.Process(target=Fitting,args = [1])
process3 = multiprocessing.Process(target=Fitting,args = [2])
process4 = multiprocessing.Process(target=Fitting,args = [3])


process1.start()
process2.start()
process3.start()
process4.start()

process1.join()
process2.join()
process3.join()
process4.join()


