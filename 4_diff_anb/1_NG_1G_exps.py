
import numpy as np
from scipy.optimize import curve_fit

from astropy.io import fits
import matplotlib.pyplot as plt
import time

t1 = time.time()

def objective(t,a,b):

    return a * (np.exp(-b*t)-1)

time1 = np.arange(0,13,1)

diff_avs = []
diff_bvs = []
count_values = []

er_a_channel_level = []
er_b_channel_level = []

for channel in range(0,4,1):
    avs = []
    bvs = []
    counts = []

    er_a=[]
    er_b=[]
    
    dc_long = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))
    dc_short = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_1G/AvgDataCube_1G_{}.npy'.format(channel))

    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(13)
    
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    count_range = 50
    for flux_bin in range(500,60000,count_range):
        pix_mask = (dc_long[-1] > (flux_bin-count_range/2)) & (dc_long[-1] < (flux_bin+count_range/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc_long[10,pix_mask])
        delta_change_cube_long = np.nanmedian((dc_long[10:10+13,pix_mask]-dc_long[9,pix_mask])*100/dc_long[9,pix_mask],axis=1)
        delta_change_cube_short = np.nanmedian((dc_short[17:,pix_mask]-dc_short[16,pix_mask])*100/dc_short[16,pix_mask],axis=1)
        print(np.shape(delta_change_cube_long),np.shape(delta_change_cube_short))
        time1 = np.arange(0,13,1)
        if not np.isnan(avg_count) and delta_change_cube_long[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
       #     print(np.shape(Sigma))
            try:
                fitting_variables_long, error1 = curve_fit(objective,time1,delta_change_cube_long,np.array([1.0,0.5]),sigma=Sigma)
                fitting_variables_short, error2 = curve_fit(objective,time1,delta_change_cube_short,np.array([1.0,0.5]),sigma=Sigma)
            
                a_long,b_long = fitting_variables_long
                a_short,b_short = fitting_variables_short

                da1 = np.sqrt(error1[0][0])
                db1 = np.sqrt(error1[1][1])
                da2 = np.sqrt(error2[0][0])
                db2 = np.sqrt(error2[1][1])

            
                diff_a = (a_long-a_short)
                diff_b = (b_long-b_short)
                print(diff_a,diff_b)
                avs.append(diff_a)
                bvs.append(diff_b)
                counts.append(avg_count)


                Da = (da1+da2)
                Db = (db1+db2)
                er_a.append(Da)
                er_b.append(Db)
            except Exception as e:
                print(e)
                pass
    t2 = time.time()
    print("channel",channel,"done",(t2-t1)//60,'min',int(t2-t1)%60,'s')
    diff_avs.append(avs)
    diff_bvs.append(bvs)
    count_values.append(counts)
    er_a_channel_level.append(er_a)
    er_b_channel_level.append(er_b)

    
color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

plt.figure(figsize=(16,9))

for chan in range(4):
    plt.errorbar(count_values[chan],diff_avs[chan],yerr=er_a_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))
    #print(len(count_values[chan]))
plt.plot([0,35000],[0,0],'--',color='black')
plt.title('change of a between long and short exposure(19/6)',fontdict={'fontsize':18})
plt.legend()

plt.xlabel('Counts(ADU)',fontsize=18)
plt.ylabel('$a_l-a_s$')
plt.ylim(-0.5,0.05)
plt.savefig('percentage_diff_a_ls.png')

plt.figure(figsize=(16,9))

for chan in range(4):
    plt.errorbar(count_values[chan],diff_bvs[chan],yerr=er_b_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))
plt.plot([0,35000],[0,0],'--',color='black')    
plt.title('change of b between long and short exposure(19/6)',fontdict={'fontsize':18})
plt.legend()
plt.ylim(-.1,0.2)
plt.xlabel('counts',fontsize=18)
plt.ylabel('$b_l-b_s$')
plt.savefig('percentage_diff_b_ls.png')

