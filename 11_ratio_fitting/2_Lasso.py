#import csv
import numpy as np
from scipy.optimize import curve_fit

from astropy.io import fits
import matplotlib.pyplot as plt
import time
from scipy.optimize import least_squares
from scipy import linalg

t1 = time.time()

def objective(t,a,b):

    return a*(np.exp(-b*t)-1)

def objective1(x,t,C,S):
    Objective = np.vstack((x[0]*(np.exp(-x[1]*t)-1)-C)/S,np.sqrt(x[2]*np.abs(x[0]))) 
    return np.vstack(Objective,np.sqrt(x[2]*np.abs(x[1])))


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

#time1 = np.arange(0,27,1)

a_values_channel_level = []
b_values_channel_level = []
count_values_channel_level = []

er_a_channel_level = []
er_b_channel_level = []



for channel in range(0,4,1):
    plt.figure()
    a_values_image_level = []
    b_values_image_level = []
    count_values_image_level = []
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
        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1)
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
         
            try:
                if len(bvs)==0:
                    b_ini = 0.5
                else:
                    b_ini = 0.5
                
                fitting_data_variables, error = curve_fit(objective,time1,delta_change_cube,np.array([1.0,0.5]),sigma=Sigma)
                a,b = fitting_data_variables
                da = np.sqrt(error[0][0])
                db = np.sqrt(error[1][1])
             
                avs.append(a)
                bvs.append(b)
                er_a.append(da)
                er_b.append(db)
                counts.append(avg_count)
                print('ch={0},S_a = {1}, S_b = {2},a={3},b={4}'.format(channel,da,db,a,b))
            except Exception as e:
                print(e)

                
                print('{0},{1}'.format(avg_count, channel,'cannot be done'))
             
    t2 = time.time()
    print("channel",channel,"done with initials",(t2-t1)//60,'min',int(t2-t1)%60,'s Curve fit')
    a_values_channel_level.append(avs)
    b_values_channel_level.append(bvs)
    count_values_channel_level.append(counts)
    er_a_channel_level.append(er_a)
    er_b_channel_level.append(er_b)



#Least Squares method

print('Curve_fit done, least square starts')

a_values_channel_level_lsq = []
b_values_channel_level_lsq = []
count_values_channel_level_lsq = []

er_a_channel_level_lsq = []
er_b_channel_level_lsq = []

for channel in range(0,4,1):
    plt.figure()
    a_values_image_level = []
    b_values_image_level = []
    count_values_image_level = []
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
        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1)
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
         
            try:
                if len(bvs)==0:
                    b_ini = 0.5
                else:
                    b_ini = 0.5
                
                fitting_data_variables= least_squares(objective1,x0=np.array([1.0,0.5,1.0]),args=(time1,delta_change_cube,Sigma))
                error = calculate_cov(fitting_data_variables)
                a,b,lambd = fitting_data_variables.x
                da = np.sqrt(error[0][0])
                db = np.sqrt(error[1][1])
                avs.append(a)
                bvs.append(b)
                er_a.append(da)
                er_b.append(db)
                counts.append(avg_count)
                print('ch={0},S_a = {1}, S_b = {2},a={3},b={4}'.format(channel,da,db,a,b),'lambda=',lambd)
            except Exception as e:
                print(e)

                
                print('{0},{1}'.format(avg_count, channel,'cannot be done'))
             
    t2 = time.time()
    print("channel",channel,"done with initials",(t2-t1)//60,'min',int(t2-t1)%60,'s for LSq')
    a_values_channel_level_lsq.append(avs)
    b_values_channel_level_lsq.append(bvs)
    count_values_channel_level_lsq.append(counts)
    er_a_channel_level_lsq.append(er_a)
    er_b_channel_level_lsq.append(er_b)

 


    




print('parametrization done')
#fig, axs = plt.subplots(4,figsize=(16,9))
plt.title('checking each exposure:a')
color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

#print(a_values_channel_level)



plt.figure(figsize=(12,6))
for chan in range(4):
     #plt.plot(count_values_channel_level[chan],a_values_channel_level[chan],'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))
     plt.errorbar(count_values_channel_level[chan],a_values_channel_level[chan],yerr=er_a_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))
     plt.errorbar(count_values_channel_level_lsq[chan],a_values_channel_level_lsq[chan],yerr=er_a_channel_level_lsq[chan],fmt='<',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='LSq Channel {}'.format(chan+1))
plt.title("a vs count: Dataset 4(19/6/22),  $a(e^{-bt}-1)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
plt.ylim(0.0,0.65)
plt.ylabel('a(%)',fontsize=18)
plt.savefig("DS5_Expo_check_a_exp_l.png")



plt.figure(figsize=(12,6))
for chan in range(4):
     
     plt.errorbar(count_values_channel_level[chan],b_values_channel_level[chan],yerr=er_b_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))
     plt.errorbar(count_values_channel_level_lsq[chan],b_values_channel_level_lsq[chan],yerr=er_b_channel_level_lsq[chan],fmt='<',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='LSq Channel {}'.format(chan+1))
     
plt.title("b vs count: Dataset 4(19/6/22), $a(e^{-bt}-1)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
plt.ylim(0.0,0.35)
plt.ylabel('b(per readout)',fontsize=18)
plt.savefig("DS5_Expo_check_b_exp_l.png")

plt.figure(figsize=(12,6))
for chan in range(4):
     #plt.plot(count_values_channel_level[chan],b_values_channel_level[chan],'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))            
     plt.errorbar(a_values_channel_level[chan],b_values_channel_level[chan],xerr=er_a_channel_level[chan],yerr=er_b_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))
     plt.errorbar(a_values_channel_level_lsq[chan],b_values_channel_level_lsq[chan],xerr=er_a_channel_level_lsq[chan],yerr=er_b_channel_level_lsq[chan],fmt='<',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='LSq Channel {}'.format(chan+1))
plt.title("b vs a: Dataset 4(19/6/22), $a(e^{-bt}-1)$",fontdict={'fontsize': 18})
plt.xlabel('a (%)',fontsize=18)
plt.legend()
plt.xlim(0.0,0.65)
plt.ylim(0.0,0.35)
plt.ylabel('b(per readout)',fontsize=18)
plt.savefig("DS5_Expo_check_b_vs_a.png")



