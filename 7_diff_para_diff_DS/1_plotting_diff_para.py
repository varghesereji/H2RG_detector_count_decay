import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time

def objective(t,a,b):
    return a*(np.exp(-b*t)-1)

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral'] 

t1 = time.time()
for channel in range(0,4,1):
    Davs = []
    Dbvs = []
    Count = []

    Der_a = []
    Der_b = []

    dc_3 = np.load('/home/varghese/Desktop/DP1/Dataset_3/long_time_expo_avg/Avg_long_time_expo_{}.npy'.format(channel))

    dc_5 = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))


    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{\
}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)

    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]

    interval = 100

    for flux_bin in range(500,60000,interval):

        pix_mask = (dc_3[-1] >(flux_bin-interval/2)) & (dc_3[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)

        avg_count = np.nanmedian(dc_3[10,pix_mask])


        delta_change_cube_1 = np.nanmedian((dc_3[10:,pix_mask]-dc_3[10,pix_mask])*100/dc_3[10,pix_mask],axis=1)
        delta_change_cube_2 = np.nanmedian((dc_5[10:,pix_mask]-dc_5[10,pix_mask])*100/dc_5[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube_1,axis=0),1)

        if not np.isnan(avg_count) and delta_change_cube_1[-1]<0 and delta_change_cube_2[-1] < 0:

            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))


            try:

   
                fitting_data_variables_1, error_1 = curve_fit(objective,time1,delta_change_cube_1,np.array([1.0,0.5]),sigma=Sigma)
                a1,b1 = fitting_data_variables_1
                da1 = np.sqrt(error_1[0][0])
                db1 = np.sqrt(error_1[1][1])

                
                fitting_data_variables_2, error_2 = curve_fit(objective,time1,delta_change_cube_2,np.array([1.0,0.5]),sigma=Sigma)

                a2,b2 = fitting_data_variables_2
                da2 = np.sqrt(error_2[0][0])
                db2 = np.sqrt(error_2[1][1])

                diff_a = (a1-a2)*100/a2
                diff_b = (b1-b2)*100/b2

                Er_a = (da1+da2)*100/a2 + da2*(a1-a2)*100/a2**2
                Er_b = (db1+db2)*100/b2 + db2*(b1-b2)*100/(b2**2)

                Davs.append(diff_a)
                Dbvs.append(diff_b)
                Count.append(avg_count)
                Der_a.append(Er_a)
                Der_b.append(Er_b)
            except Exception as e:
                print(e)

    t2 = time.time()
    np.save('Diff_a_{}.npy'.format(channel),Davs)
    np.save('Diff_b_{}.npy'.format(channel),Dbvs)         
    np.save('Count_{}.npy'.format(channel),Count)
    np.save('Er_a_{}.npy'.format(channel),Der_a)
    np.save('Er_b_{}.npy'.format(channel),Der_b)


    print("channel",channel,"done with initials",(t2-t1)//60,'min',int(t2-t1)%60,'s')

                
                
