import csv
import numpy as np
from scipy.optimize import curve_fit

from astropy.io import fits
import matplotlib.pyplot as plt
import time

t1 = time.time()

def P1(x):
    return x
def P2(x):
    return (1/2)*(3*t**2-1)

def objective(t,c,d):
    
    print(t)
    return -c*P1(t)+d*(2*P2(t)/3+1/3)

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
    fields = ['cro','avg_count','a','b','fit']
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    interval = 100
    for flux_bin in range(500,60000,interval):
        #rows = []
        
        cro+=1
        pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc[10,pix_mask])
        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / np.size(delta_change_cube,axis=0)
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
            #Sigma = Sigma*delta_change_cube/(np.nanmedian(dc[3:,pix_mask]-dc[3,pix_mask])*np.sqrt(np.sum(pix_mask)))
            #print(Sigma,channel,np.sum(pix_mask))
            sub_rows = []
            try:
                if len(bvs)==0:
                    b_ini = 0.5
                else:
                    b_ini = 0.5
                
                fitting_data_variables, error = curve_fit(objective,time1,delta_change_cube,np.array([1.0,0.5]),sigma=Sigma)
                a,b = fitting_data_variables
                da = np.sqrt(error[0][0])
                db = np.sqrt(error[1][1])
                #df = np.sqrt(error[2][2])
                avs.append(a)
                bvs.append(b)
                er_a.append(da)
                er_b.append(db)
                counts.append(avg_count)
               # print(Sigma)
                #print(channel,cro-1,'avg count=',avg_count,'# pixels=',np.sqrt(np.sum(pix_mask)),'Read error=',np.median(np.sqrt(Sigma2.flatten())),'Fin Sigma=',np.median(np.sqrt(Sigma2.flatten()))*100/(avg_count*np.sqrt(np.sum(pix_mask))),'fit with initial')
                print('ch={0},S_a = {1}, S_b = {2},a={3},b={4}'.format(channel,da,db,a,b))
                sub_rows = [cro-1,avg_count,a,b,'yes']
                rows.append(sub_rows)
            except Exception as e:
                print(e)
                #if avg_count != nan:
                
                print('{0},{1}'.format(avg_count, channel,'cannot be done'))

                sub_rows = [cro-1,avg_count,e,'nil','no']
                rows.append(sub_rows)
    t2 = time.time()
    print("channel",channel,"done with initials",(t2-t1)//60,'min',int(t2-t1)%60,'s')
    #a_values_image_level.append(avs)
    #b_values_image_level.append(bvs)
    #count_values_image_level.append(counts)
    a_values_channel_level.append(avs)
    b_values_channel_level.append(bvs)
    count_values_channel_level.append(counts)
    er_a_channel_level.append(er_a)
    er_b_channel_level.append(er_b)
    #plt.legend()
    #plt.title('Wrong fittings')
    #plt.savefig('wrong_curves_ch{}.png'.format(channel))
    #print(rows)
    filename = "Poly_range100_Channel_{}_data.csv".format(channel)
    with open(filename,'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)
    print("CSV file",filename,"Saved")
print('parametrization done')
#fig, axs = plt.subplots(4,figsize=(16,9))
plt.title('checking each exposure:a')
color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

#print(a_values_channel_level)

plt.figure(figsize=(12,6))
for chan in range(4):
     #plt.plot(count_values_channel_level[chan],a_values_channel_level[chan],'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))
     plt.errorbar(count_values_channel_level[chan],a_values_channel_level[chan],yerr=er_a_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))#,'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))
plt.title("a vs count: Dataset 4(19/6/22),  $a(e^{-bt}-1)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
plt.ylim(0.0,0.65)
plt.ylabel('a(%)',fontsize=18)
plt.savefig("DS5_Expo_check_a_exp_l.png")



plt.figure(figsize=(12,6))
for chan in range(4):
     #plt.plot(count_values_channel_level[chan],b_values_channel_level[chan],'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))
     plt.errorbar(count_values_channel_level[chan],b_values_channel_level[chan],yerr=er_b_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))#,'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))
plt.title("b vs count: Dataset 4(19/6/22), $a(e^{-bt}-1)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
plt.ylim(0.0,0.35)
plt.ylabel('b(per readout)',fontsize=18)
plt.savefig("DS5_Expo_check_b_exp_l.png")

plt.figure(figsize=(12,6))
for chan in range(4):
     #plt.plot(count_values_channel_level[chan],b_values_channel_level[chan],'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))            
     plt.errorbar(a_values_channel_level[chan],b_values_channel_level[chan],xerr=er_a_channel_level[chan],yerr=er_b_channel_level[chan],fmt='.',color=color_list[chan],ecolor=color_list[4+chan],elinewidth=1,capsize=10,label='Channel {}'.format(chan+1))#,'.',color=color_list[chan],alpha=0.7,label='{}'.format(chan))           
plt.title("b vs a: Dataset 4(19/6/22), $a(e^{-bt}-1)$",fontdict={'fontsize': 18})
plt.xlabel('a (%)',fontsize=18)
plt.legend()
plt.xlim(0.0,0.65)
plt.ylim(0.0,0.35)
plt.ylabel('b(per readout)',fontsize=18)
plt.savefig("DS5_Expo_check_b_vs_a.png")



