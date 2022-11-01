#import csv
import numpy as np
from scipy.optimize import curve_fit
import multiprocessing
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from scipy.optimize import least_squares
from scipy import linalg

t1 = time.time()

def objective(t,b):

    return (np.exp(-b*t))

def objective1(x,t,C,S,lamb):
    Objective = np.append(((np.exp(-x[0]*t)-1)-C)/S,np.sqrt(lamb*np.abs(x[0])))
    Objective = np.append(Objective,np.sqrt(lamb*np.abs(x[0])))

    return Objective 


def calculate_cov(scipylsq_res,absolute_sigma=True):
    _,s,VT=linalg.svd(scipylsq_res.jac,full_matrices=False)
    threshold = np.finfo(float).eps * max(scipylsq_res.jac.shape)
    s=s[s>threshold]
    VT = VT[:s.size]
    pcov = np.dot(VT.T/s**2,VT)
    if not absolute_sigma:
        ysize = len(scipylsq_res.fun)
        cost = 2*scipylsq_res.cost
        popt = scipylsq_res.x
        if ysize > len(popt):
            s_sq = cost / (ysize - len(popt))
            pcov = pcov*s_sq
    return pcov



def Fitting(lamb,channel):

    avs = []
    bvs = []

    er_a = []
    er_b = []
    counts = []
    cro = 0
    dc = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))
    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)
    # print(Sigma)
    rows = []

    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    interval = 100
    for flux_bin in range(500,60000,interval):
        #rows = []

        cro+=1
        pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc[10,pix_mask])
        delta_change_cube = np.nanmedian((dc[10:,pix_mask])/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1)
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

            try:
                if len(bvs)==0:
                    b_ini = 0.5
                else:
                    b_ini = 0.5

                fitting_data_variables= least_squares(objective1,x0=np.array([0.5]),args=(time1,delta_change_cube,Sigma,lamb))
                error = calculate_cov(fitting_data_variables)
                b = fitting_data_variables.x
                #da = np.sqrt(error[0][0])
                db = np.sqrt(error[0][0])
                print("b=",b)
                #avs.append(a)
                bvs.append(b)
                #er_a.append(da)
                er_b.append(db)
                counts.append(avg_count)


                plt.figure(figsize=(16,9))
                plt.plot(time1,delta_change_cube,'o-',color='blue')
                plt.plot(time1,objective(time1,b))
                plt.xlabel('time')
                plt.ylabel('count difference')
                plt.title('b = {1}'.format(b))
                plt.savefig('Fitting/Plot_{0}_{1}.png'.format(channel,flux_bin))
                
             #   print('ch={0},S_a = {1}, S_b = {2},a={3},b={4}'.format(channel,da,db,a,b),'lamb=',lamb)
            except Exception as e:
                print(e)


                print('{0},{1}'.format(avg_count, channel,'cannot be done'))

    return bvs,er_b,counts

def saving_files(lamb,channel):
    bvs,er_b,counts = Fitting(lamb,channel)

    #np.save('a_values/avs_{0}.npy'.format(channel),np.array(avs))
    np.save('b_values/bvs_{0}.npy'.format(channel),np.array(bvs))
    np.save('Count_values/avs_{0}.npy'.format(channel),np.array(counts))
    #np.save('era_values/eravs_{0}.npy'.format(channel),np.array(er_a))
    np.save('erb_values/erbvs_{0}.npy'.format(channel),np.array(er_b))
       
    t2 = time.time()
    print("lambda=",lamb,"channel",channel,"done with initials",(t2-t1)//60,'min',int(t2-t1)%60,'s for LSq')
    



#for lamb in range(10000000,110000000,1000000):
    
lamb = 0

process1 = multiprocessing.Process(target=saving_files,args = [lamb,0])
process2 = multiprocessing.Process(target=saving_files,args = [lamb,1])
process3 = multiprocessing.Process(target=saving_files,args = [lamb,2])
process4 = multiprocessing.Process(target=saving_files,args = [lamb,3])

process1.start()
process2.start()
process3.start()
process4.start()

process1.join()
process2.join()
process3.join()
process4.join()


 


    
 


