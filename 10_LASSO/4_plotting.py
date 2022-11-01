import numpy as np
import matplotlib.pyplot as plt

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']


    

plt.figure(figsize=(16,9))
for channel in range(0,4,1):

    avs = np.load('a_values/avs_{0}.npy'.format(channel))
    bvs = np.load('b_values/bvs_{0}.npy'.format(channel))
    #coun = np.load('/home/varghese/Desktop/Least_sq/a_values/avs_{0}_{1}.npy'.format(channel,lamb))
    er_a = np.load('era_values/eravs_{0}.npy'.format(channel))
    er_b = np.load('erb_values/erbvs_{0}.npy'.format(channel))

    plt.errorbar(avs,bvs,xerr=er_a,yerr=er_b,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
    print('channel = {} done'.format(channel))
plt.xlabel('a(%)')
plt.ylabel('b(per readouts)')
plt.title('Correlation btw a and b')


plt.savefig('Plot_correlation.png')

    
        

    
