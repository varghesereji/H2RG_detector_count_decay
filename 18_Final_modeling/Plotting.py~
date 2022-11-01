import numpy as np
import matplotlib.pyplot as plt


color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

plt.figure()
for channel in range(4):
    counts = np.load('new_count_{}.npy'.format(channel))
    pvs = np.load('new_pvs_{}.npy'.format(channel))
    err = np.load('new_p_err_{}.npy'.format(channel))
    print(np.median(err))
    plt.errorbar(counts,pvs,yerr=err,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
plt.ylim(-0.05,0.6)
plt.xlabel('counts')
plt.ylabel('p')
plt.title('p vs counts from exponential')
plt.savefig('1_p_from_exp.png')
    
