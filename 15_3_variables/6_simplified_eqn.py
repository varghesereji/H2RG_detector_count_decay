

import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import multiprocessing

t1 = time.time()

def P1(x):
    return x

def P2(x):
    return (1/2)*(3*x**2-1)


def objective(t,c,d):
#    print(t)
    p1 = P1(t)
    p2 = P2(t)

    return c*(t+1)+d*(t**2-1)


def Fitting(channel):
    cvs = []
    dvs = []
    er_c = []
    er_d = []
    counts = []

    dc = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))

    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten()))+np.zeros(20)
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:512*(channel+1)]

    interval = 100
    for flux_bin in range(500,60000,interval):

        pix_mask = (dc[-1] > (flux_bin-interval//2)) & (dc[-1] < (flux_bin+interval//2)) & (fiber_mask == 1)

        avg_count = np.nanmedian(dc[10,pix_mask])

        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / (np.size(delta_change_cube,axis=0)/2)-1

        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

            try:
                fitting_data_variables, error = curve_fit(objective,time1,delta_change_cube,np.array([1.0,0.5]),sigma = Sigma)
                print(error)
                c,d = fitting_data_variables

                dc1 = np.sqrt(error[0][0])
                dd = np.sqrt(error[1][1])

                plt.figure(figsize=(16,9))
                plt.plot(time1,delta_change_cube,'o-',color='blue')
                plt.plot(time1,objective(time1,c,d),'--',color='green')
                plt.xlabel('time')
                plt.ylabel('count')
                plt.title('Count level:{}'.format(flux_bin))
                plt.savefig('Ch_{0}/Fitted_{1}.png'.format(channel,flux_bin))
                
                cvs.append(c)
                dvs.append(d)
                counts.append(avg_count)
                er_c.append(dc1)
                er_d.append(dd)
            except Exception as e:
                print(e)

                print('{0},{1}'.format(avg_count,channel),'cannot be done')

    t2 = time.time()

    print("channel", channel, 'done with initials',(t2-t1)//60, 'min', int(t2-t1)%60,'s')

    np.save('cvs/cvs_{}.npy'.format(channel),cvs)
    np.save('dvs/dvs_{}.npy'.format(channel),dvs)
    np.save('counts/counts_{}.npy'.format(channel),counts)
    np.save('er_c/er_c_{}.npy'.format(channel),er_c)
    np.save('er_d/er_d_{}.npy'.format(channel),er_d)



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




#for channel in range(0,4,1):
#    Fitting(channel)
