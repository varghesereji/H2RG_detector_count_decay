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
    dc_file_0G = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(q)
    dc_file_1G = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_1G/AvgDataCube_1G_{}.npy'.format(q)
    dc_file_2G = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_2G/AvgDataCube_2G_{}.npy'.format(q)
    dc_0G= np.load(dc_file_0G)
    dc_1G=np.load(dc_file_1G)
    dc_2G = np.load(dc_file_2G)

    print(q,"started")
    Dist_slope1 = None
    Dist_slope2 = None

    for flux_bin in range(1000,50000,50):
        
        pix_mask = (dc_0G[-1] > (flux_bin - 25)) & (dc_0G[-1] < (flux_bin + 25))
        #avg_count = np.nanmedian(dc[10,pix_mask])

        tr_0G = np.nanmedian(dc_0G[3:8,pix_mask],axis=1)
        tr_2G = np.nanmedian(dc_2G[17:22,pix_mask],axis=1)
        tr_1G = np.nanmedian(dc_1G[10:15,pix_mask],axis=1)
        len_tr = np.size(dc_0G[2:7,pix_mask],axis=1)
        #print(len_tr)
        d_0G = np.diff(tr_0G,axis=0)
        d_1G = np.diff(tr_1G,axis=0)
        d_2G = np.diff(tr_2G,axis=0)
        
        slope1 = d_0G / d_1G
        slope2 = d_0G / d_2G

        if Dist_slope1 is None:
            Dist_slope1 = slope1
            Dist_slope2 = slope2
        else:
            Dist_slope1 = np.vstack((Dist_slope1,slope1))
            Dist_slope2 = np.vstack((Dist_slope2,slope2))
        #plt.plot(slope,color=color_list[q],alpha=0.3)
    Dist_slope1 = np.transpose(Dist_slope1)
    Dist_slope2=np.transpose(Dist_slope2)
    t=0
    print('chan','\t time','\t mu1','\t mu2','\t SD1','\t SD2','\t Dev1','\t Dev2')
    for slopes in range(len(Dist_slope1)):
        mu1=np.nanmean(Dist_slope1[slopes])
        mu2 = np.nanmean(Dist_slope2[slopes])
        num_of_points = np.sum(np.isfinite(Dist_slope1[slopes]))
        SD1 = np.nanstd(Dist_slope1[slopes])/np.sqrt(num_of_points)
        SD2 = np.nanstd(Dist_slope2[slopes])/np.sqrt(num_of_points)
        #med = np.nanmedian(slopes2)
        diff1 = np.abs(1-mu1)
        diff2 = np.abs(1-mu2)
        n_sigma1=int(diff1/SD1)+1
        n_sigma2 = int(diff2 / (SD2)) +1 
        print(q,t,mu1,mu2,SD1,SD2, n_sigma1, n_sigma2,)
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
