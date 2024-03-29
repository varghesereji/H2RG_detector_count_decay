

import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import multiprocessing

t1 = time.time()

def v1(x):
    return x

def v2(x):
    return ((-3/8)*x+x**2/2)


def objective(t,p,q):
#    print(t)
    V1 = v1(t)
    V2 = v2(t)

    return -p*V1+q*V2

def objective1(t,a,b):
    return a*(np.exp(-b*t)-1)


def Fitting(channel):
    pvs = []
    qvs = []
    er_p = []
    er_q = []
    counts = []

    avs = []
    bvs = []
    er_a = []
    er_b = []
    
    dc = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))

    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten()))+np.zeros(20)
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:512*(channel+1)]

    interval = 250
    for flux_bin in range(500,60000,interval):

        pix_mask = (dc[-1] > (flux_bin-interval//2)) & (dc[-1] < (flux_bin+interval//2)) & (fiber_mask == 1)

        avg_count = np.nanmedian(dc[10,pix_mask])

        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / (np.size(delta_change_cube,axis=0))

        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

            try:
                fitting_data_variables, error = curve_fit(objective,time1,delta_change_cube,np.array([1.0,0.5]),sigma = Sigma)

                p,q = fitting_data_variables

                dp = np.sqrt(error[0][0])
                dq = np.sqrt(error[1][1])

                fitting_ab,error1 = curve_fit(objective1,time1,delta_change_cube,np.array([0.01,0.01]),sigma=Sigma)
                a,b = fitting_ab
                da = np.sqrt(error1[0][0])
                db = np.sqrt(error1[1][1])
                np.save('Ch_{0}/Decay_curve_{1}_{2}.npy'.format(channel,channel,flux_bin),delta_change_cube)

                pvs.append(p)
                qvs.append(q)
                counts.append(flux_bin)
                er_p.append(dp)
                er_q.append(dq)
                avs.append(a)
                bvs.append(b)
                er_a.append(da)
                er_b.append(db)
            except Exception as e:
                print(e)

                print('{0},{1}'.format(avg_count,channel),'cannot be done')
 
    t2 = time.time()

    print("channel", channel, 'done with initials',(t2-t1)//60, 'min', int(t2-t1)%60,'s')

    np.save('pvs/pvs_{}.npy'.format(channel),pvs)
    np.save('qvs/qvs_{}.npy'.format(channel),qvs)
    np.save('counts/counts_{}.npy'.format(channel),counts)
    np.save('er_p/er_p_{}.npy'.format(channel),er_p)
    np.save('er_q/er_q_{}.npy'.format(channel),er_q)
    np.save('avs/avs_{}.npy'.format(channel),avs)
    np.save('bvs/bvs_{}.npy'.format(channel),bvs)
    np.save('er_a/er_a_{}.npy'.format(channel),er_a)
    np.save('er_b/er_b_{}.npy'.format(channel),er_b)
    


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
