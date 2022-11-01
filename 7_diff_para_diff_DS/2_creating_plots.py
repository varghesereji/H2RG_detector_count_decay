
import numpy as np
import matplotlib.pyplot as plt


color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']



plt.figure(figsize=(16,9))
for channel in range(4):
    a_vs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/7_diff_para_diff_DS/Diff_a_{}.npy'.format(channel))

    count = np.load('/home/varghese/Desktop/DS5_DP1/Codes/7_diff_para_diff_DS/Count_{}.npy'.format(channel))

    Err_a = np.load('/home/varghese/Desktop/DS5_DP1/Codes/7_diff_para_diff_DS/Er_a_{}.npy'.format(channel))

    plt.errorbar(count,a_vs,yerr=Err_a,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
plt.ylim(-50,50)
plt.title('Difference of $a$ for 13/5 Dataset and 19/6 Dataset')
plt.xlabel('Counts(ADU)')
plt.ylabel('$\dfrac{ a_{13/5}-a_{19/6} }{a_{19/6}}X 100$')
plt.legend()
plt.grid()
plt.savefig('Diff_DS_a.png')
    

plt.figure(figsize=(16,9))
for channel in range(4):
    b_vs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/7_diff_para_diff_DS/Diff_b_{}.npy'.format(channel))

    count = np.load('/home/varghese/Desktop/DS5_DP1/Codes/7_diff_para_diff_DS/Count_{}.npy'.format(channel))

    Err_b = np.load('/home/varghese/Desktop/DS5_DP1/Codes/7_diff_para_diff_DS/Er_b_{}.npy'.format(channel))

    plt.errorbar(count,b_vs,yerr=Err_b,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
plt.ylim(-50,50)
plt.xlabel('Counts(ADU)')
plt.ylabel('$\dfrac{b_{13/5}-b_{19/6}}{b_{19/6}} X 100$(per readout)')
plt.title('Difference of $b$ for 13/5 Dataset and 19/6 Dataset')
plt.legend()
plt.grid()
plt.savefig('Diff_DS_b.png')
    
