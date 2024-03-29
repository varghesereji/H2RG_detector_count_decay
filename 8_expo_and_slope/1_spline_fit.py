#\65;6800;1c#This script is to compare the slopes of lobes of dataset on 19/6 using spline fit


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import UnivariateSpline

color_list = ['k','b','g','r']

def objective(t,a,b):
    return a*(np.exp(-b*t)-1)/100


a_min = 0.08
a_max = 0.25
b_max = 0.05
b_min = 0.1


plt.figure(figsize=(16,9))
for q in range(4):
    dc_0G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(q))
    dc_1G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_1G/AvgDataCube_1G_{}.npy'.format(q))
    dc_2G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_2G/AvgDataCube_2G_{}.npy'.format(q))

   # Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*q:(512*(q+1))]

    print(q,'started')
    
    Dist_slope_1=[]
    Dist_slope_2=[]

    selected_flux_list_1 = []
    selected_flux_list_2 = []

    count_range = 100
    
    for flux_bin in range(500,60000,count_range):
        pix_mask = (dc_0G[-1] > (flux_bin-count_range/2)) & (dc_0G[-1] < (flux_bin+count_range/2)) & (fiber_mask == 1)

        tr_0G = np.nanmedian(dc_0G[3:8,pix_mask],axis=1)
        tr_1G = np.nanmedian(dc_1G[10:15,pix_mask],axis=1)
        tr_2G = np.nanmedian(dc_2G[17:22,pix_mask],axis=1)

        try:
            tr_0G_spl = UnivariateSpline(tr_0G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
            tr_1G_spl = UnivariateSpline(tr_1G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
            tr_2G_spl = UnivariateSpline(tr_2G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
#		fitting_data_variables,error = curve_fit(objective,time1,delta_change_cube,np.array([1.0,0.5]),sigma=Sigma)
 #               a,b = fitting_data_variables

  #              Gc = 1/(objective(T1,a,b)+1)
            d_0G = 1/tr_0G_spl(tr_0G,nu=1)
            d_1G = 1/tr_1G_spl(tr_1G,nu=1)
            d_2G = 1/tr_2G_spl(tr_2G,nu=1)

            slope_1 = d_0G / d_1G
            slope_2 = d_0G / d_2G
            Dist_slope_1.append(slope_1)
            Dist_slope_2.append(slope_2)
            print(np.nanmedian(slope_1),np.nanmedian(slope_2))
#            Gain_change_expo.append(Gc)
            selected_flux_list_1.append(flux_bin)
        except Exception as e:
             print(e)

        
        
    plt.plot(selected_flux_list_1,np.median(Dist_slope_1,axis=1),'o',color=color_list[q],label='Channel {} 0G/1G'.format(q+1))
    plt.plot(selected_flux_list_1,np.median(Dist_slope_2,axis=1),'<',color=color_list[q],label='Channel {} 0G/2G'.format(q+1))


T1 = 15
T2 = 22

GC1_min = objective(T1,a_min,b_min)+1
GC1_max = objective(T1,a_max,b_max)+1
GC2_min = objective(T2,a_min,b_min)+1
GC2_max = objective(T2,a_max,b_max)+1

#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC1_min,'--',color='k',label='Time T1(End of 1G) a={0}, b={1}'.format(a_min,b_min))
#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC1_max,':',color='k',label='Time T1(End of 1G) a={0}, b={1}'.format(a_max,b_max))
#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC2_min,'--',color='b',label='Time T2(End of 2G) a = {0}, b={1}'.format(a_min,b_min))
#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC2_max,':',color='b',label='Time T1(End of 1G) a={0}, b={1}'.format(a_max,b_max))

plt.ylim(0.9975,1.0025)
plt.grid()
plt.legend()
plt.xlabel('Count',fontsize=18)
plt.ylabel('Ratio of slopes',fontsize=18)
plt.title('Ratio of slopes by spline fit',fontsize=18)
plt.savefig('Splive_fit_slope_vs_count.png')
    









plt.figure(figsize=(16,9))
for q in range(4):
    dc_0G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(q))
    dc_1G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_1G/AvgDataCube_1G_{}.npy'.format(q))
    dc_2G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_2G/AvgDataCube_2G_{}.npy'.format(q))

   # Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*q:(512*(q+1))]

    print(q,'started')
    
    Dist_slope_1=[]
    Dist_slope_2=[]

    selected_flux_list_1 = []
    selected_flux_list_2 = []

    count_range = 100
    
    for flux_bin in range(500,60000,count_range):
        pix_mask = (dc_0G[-1] > (flux_bin-count_range/2)) & (dc_0G[-1] < (flux_bin+count_range/2)) & (fiber_mask == 1)

        tr_0G = np.nanmedian(dc_0G[3:8,pix_mask],axis=1)
        tr_1G = np.nanmedian(dc_1G[10:15,pix_mask],axis=1)
        tr_2G = np.nanmedian(dc_2G[17:22,pix_mask],axis=1)

        try:
            tr_0G_spl = UnivariateSpline(tr_0G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
            tr_1G_spl = UnivariateSpline(tr_1G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
            tr_2G_spl = UnivariateSpline(tr_2G,np.arange(len(tr_0G)),k=2,s=0,ext=0)
#		fitting_data_variables,error = curve_fit(objective,time1,delta_change_cube,np.array([1.0,0.5]),sigma=Sigma)
 #               a,b = fitting_data_variables

  #              Gc = 1/(objective(T1,a,b)+1)
            d_0G = 1/tr_0G_spl(tr_0G,nu=1)
            d_1G = 1/tr_1G_spl(tr_1G,nu=1)
            d_2G = 1/tr_2G_spl(tr_2G,nu=1)

            slope_1 = d_0G / d_1G
            slope_2 = d_0G / d_2G
            Dist_slope_1.append(slope_1)
            Dist_slope_2.append(slope_2)
            print(np.nanmedian(slope_1),np.nanmedian(slope_2))
#            Gain_change_expo.append(Gc)
            selected_flux_list_1.append(flux_bin)
        except Exception as e:
             print(e)

        
        
    plt.plot(np.median(Dist_slope_1,axis=1),np.median(Dist_slope_1,axis=1),'o',color=color_list[q],label='Channel {} 0G/1G'.format(q+1))
    #plt.plot(selected_flux_list_1,np.median(Dist_slope_2,axis=1),'<',color=color_list[q],label='Channel {} 0G/2G'.format(q+1))


T1 = 15
T2 = 22

GC1_min = objective(T1,a_min,b_min)+1
GC1_max = objective(T1,a_max,b_max)+1
GC2_min = objective(T2,a_min,b_min)+1
GC2_max = objective(T2,a_max,b_max)+1

#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC1_min,'--',color='k',label='Time T1(End of 1G) a={0}, b={1}'.format(a_min,b_min))
#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC1_max,':',color='k',label='Time T1(End of 1G) a={0}, b={1}'.format(a_max,b_max))
#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC2_min,'--',color='b',label='Time T2(End of 2G) a = {0}, b={1}'.format(a_min,b_min))
#plt.plot(selected_flux_list_1,np.zeros(np.size(selected_flux_list_1))+1/GC2_max,':',color='b',label='Time T1(End of 1G) a={0}, b={1}'.format(a_max,b_max))

plt.ylim(0.9975,1.0025)
plt.grid()
plt.legend()
plt.xlabel('Ratio of Slope 1 (NG/1G)',fontsize=18)
plt.ylabel('Ratio of Slope 2(NG/2G)',fontsize=18)
plt.title('Slope1 vs Slope2',fontsize=18)
plt.savefig('Splive_fit_slope1_vs_slope2.png')
    
