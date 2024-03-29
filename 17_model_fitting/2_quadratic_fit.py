

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time
import multiprocessing

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

t1 = time.time()

def objective(x,l,m,n):
    return l*x**2+m*x+n


def Fitting(channel):
    pvs=np.load('pvs/pvs_{}.npy'.format(channel))
    qvs=np.load('qvs/qvs_{}.npy'.format(channel))
    er_p=np.load('er_p/er_p_{}.npy'.format(channel))
    er_q=np.load('er_q/er_q_{}.npy'.format(channel))


    fitted_values,error = curve_fit(objective,pvs,qvs,np.array([0.5,0.5,0.0]),sigma=er_q)
    
    l,m,n = fitted_values
    er_l = np.sqrt(error[0][0])
    er_m = np.sqrt(error[1][1])
    er_n=np.sqrt(error[2][2])
    print(channel,l,m,n,error)

    plt.figure()
    plt.errorbar(pvs,qvs,xerr=er_p,yerr=er_q,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
    plt.plot(pvs,objective(pvs,l,m,n),'--',color=color_list[channel])
    plt.savefig("Fitted_model_{}.png".format(channel))
    print("Saved")

    
    

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


