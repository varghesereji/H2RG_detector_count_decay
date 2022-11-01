#This script is to compare the slopes of lobes of dataset on 19/6 using spline fit


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import UnivariateSpline

color_list = ['k','b','g','r']

plt.figure(figsize=(16,9))
for q in range(4):
    dc_0G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(q))
    dc_1G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_1G/AvgDataCube_1G_{}.npy'.format(q))
    dc_2G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_2G/AvgDataCube_2G_{}.npy'.format(q))

    print(q,'started')

    Dist_slope_1=None
    Dist_slope_2=None

    selected_flux_list_1 = []
    selected_flux_list_2 = []

    count_range = 100

    for flux_bin in range(1000,50000,count_range):
        pix_mask = (dc_0G[-1] > (flux_bin-count_range/2)) & (dc_0G[-1] < (flux_bin+count_range/2))

        tr_0G = np.nanmedian(dc_0G[3:8,pix_mask],axis=1)
        tr_1G = np.nanmedian(dc_1G[10:15,pix_mask],axis=1)
        tr_2G = np.nanmedian(dc_2G[17:22,pix_mask],axis=1)

        try:
            tr_0G_spl = UnivariateSpline(tr_0G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
            tr_1G_spl = UnivariateSpline(tr_1G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
            tr_2G_spl = UnivariateSpline(tr_2G,np.arange(len(tr_0G)),k=2,s=0,ext=0)

        except ValueError:
            continue
        else:
            selected_flux_list_1.append(flux_bin)

        d_0G = 1/tr_0G_spl(tr_0G,nu=1)
        d_1G = 1/tr_1G_spl(tr_1G,nu=1)
        d_2G = 1/tr_2G_spl(tr_2G,nu=1)
        
        slope_1 = d_0G / d_1G
        slope_2 = d_0G / d_2G

        if Dist_slope_1 is None:
            Dist_slope_1 = slope_1
            Dist_slope_2 = slope_2
        else:
            Dist_slope_1 = np.vstack((Dist_slope_1,slope_1))
            Dist_slope_2 = np.vstack((Dist_slope_2,slope_1))
    plt.plot(selected_flux_list_1,np.median(Dist_slope_1,axis=1),'o',color=color_list[q],label='Channel {} 0G/1G'.format(q+1))
    plt.plot(selected_flux_list_1,np.median(Dist_slope_2,axis=1),'<',color=color_list[q],label='Channel {} 0G/2G'.format(q+1))

plt.legend()
plt.xlabel('Count',fontsize=18)
plt.ylabel('Ratio of slopes',fontsize=18)
plt.title('Ratio of slopes by spline fit',fontsize=18)
plt.savefig('Splive_fit_slope_vs_count.png')
    
