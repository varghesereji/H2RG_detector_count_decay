# This script is to compare the slopes of lobes of the dataset


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

datacube_list = ['/home/varghese/Desktop/DP1/Dataset_3/long_time_expo_avg/Avg_long_time_expo_{}.npy'.format(i) for i in range(4)]

color_list = ['k','b','g','r']

#fig = plt.figure(figsize=(25,15))
#gs = fig.add_gridspec(4,4,hspace=0,wspace=0)
#axs = gs.subplots(sharex='col',sharey='row')
fig, axs = plt.subplots(4,5)

for q in range(4):
    dc_file_long = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(q)
    dc_file_short = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_2G/AvgDataCube_2G_{}.npy'.format(q)
    dc_long= np.load(dc_file_long)
    dc_short = np.load(dc_file_short)

    print(q,"started")

    Dist_slope = None

    for flux_bin in range(1000,50000,50):
        
        pix_mask = (dc_long[-1] > (flux_bin - 25)) & (dc_long[-1] < (flux_bin + 25))
        #avg_count = np.nanmedian(dc[10,pix_mask])

        tr_long = np.nanmedian(dc_long[3:8,pix_mask],axis=1)
        tr_short = np.nanmedian(dc_short[17:22,pix_mask],axis=1)
        len_tr = np.size(dc_long[2:7,pix_mask],axis=1)
        #print(len_tr)
        d_long = np.diff(tr_long,axis=0)
        d_short = np.diff(tr_short,axis=0)

        slope = d_long / d_short

        if Dist_slope is None:
            Dist_slope = slope
        else:
            Dist_slope = np.vstack((Dist_slope,slope))

        #plt.plot(slope,color=color_list[q],alpha=0.3)
    Dist_slope = np.transpose(Dist_slope)
    t=0
    for slopes in Dist_slope:
        mu = np.nanmean(slopes)
        num_of_points = np.sum(np.isfinite(slopes))
        SD = np.nanstd(slopes)/np.sqrt(num_of_points)
        med = np.nanmedian(slopes)
        diff = np.abs(1-mu)
        n_sigma = int(diff / (SD)) +1 
        print(q,t,"mu=",mu,"mu_sigma=",SD,"median=",med,'num',num_of_points, "This is in", n_sigma,"sigma")
        #axs[q,t].set_title(q)
        #axs[q,t].plot(np.sort(slopes),stats.norm.pdf(np.sort(slopes),mu,SD),color=color_list[q])
        #axs[q,t].hist(slopes,bins=25,density=True, alpha=0.6,color=color_list[q])
        #axs[q,t].plot([1,1],[0,200],'--',color='black')
        #axs[q,t].plot([mu-SD/2,mu-SD],[0,200],'--',color='red')
        #axs[q,t].plot([mu+SD/2,mu+SD],[0,200],'--',color='red')
        #axs[q,t].plot([mu,mu],[0,200],'--',color='blue')
        t=t+1

    print(q,"done")

#plt.savefig("Slope_histogram.png")
#plt.show(block=False)
