import numpy as np
import matplotlib.pyplot as plt

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']
plt.figure(figsize = (16,9))
for channel in range(4):
    
    cvs = np.load('cvs/cvs_{}.npy'.format(channel))
    dvs = np.load('dvs/dvs_{}.npy'.format(channel))
    counts = np.load('counts/counts_{}.npy'.format(channel))

    er_c = np.load('er_c/er_c_{}.npy'.format(channel))
    er_d = np.load('er_d/er_d_{}.npy'.format(channel))

    plt.errorbar(counts,cvs,yerr=er_c,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("c vs count: Dataset 4(19/6/22),  $-cP_1(t)+(d)(P_2(t)+2)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('c',fontsize=18)
plt.savefig("c_legendre.png")




plt.figure(figsize = (16,9))
for channel in range(4):
    
    cvs = np.load('cvs/cvs_{}.npy'.format(channel))
    dvs = np.load('dvs/dvs_{}.npy'.format(channel))
    counts = np.load('counts/counts_{}.npy'.format(channel))

    er_c = np.load('er_c/er_c_{}.npy'.format(channel))
    er_d = np.load('er_d/er_d_{}.npy'.format(channel))

    plt.errorbar(counts,dvs,yerr=er_d,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("d vs count: Dataset 4(19/6/22),  $-cP_1(t)+(d)(P_2(t)+2)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('d',fontsize=18)
plt.savefig("d_legendre.png")





plt.figure(figsize = (16,9))
for channel in range(4):
    
    #cvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/cvs/cvs_{}.npy'.format(channel))
    avs = np.load('avs/avs_{}.npy'.format(channel))
    counts = np.load('counts/counts_{}.npy'.format(channel))

    #er_c = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/er_c/er_c_{}.npy'.format(channel))
    er_a = np.load('er_a/er_a_{}.npy'.format(channel))

    plt.errorbar(counts,avs,yerr=er_a,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("a vs count: Dataset 4(19/6/22),  $aP_0+cP_1+dP_2$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('a',fontsize=18)
plt.savefig("a_legendre.png")




 
plt.figure(figsize = (16,9))
for channel in range(4):
    
    cvs = np.load('cvs/cvs_{}.npy'.format(channel))
    dvs = np.load('dvs/dvs_{}.npy'.format(channel))
    avs = np.load('avs/avs_{}.npy'.format(channel))
    counts = np.load('counts/counts_{}.npy'.format(channel))

    er_c = np.load('er_c/er_c_{}.npy'.format(channel))
    er_d = np.load('er_d/er_d_{}.npy'.format(channel))
    er_a = np.load('er_a/er_a_{}.npy'.format(channel))

    plt.errorbar(cvs,dvs,xerr=er_c,yerr=er_d,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("c vs d: Dataset 4(19/6/22),  $-cP_1(t)+(d)(P_2(t)+2)$",fontdict={'fontsize': 18})
plt.xlabel('c',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('d',fontsize=18)
plt.savefig("c_d.png")


 
plt.figure(figsize = (16,9))
for channel in range(4):
    
    cvs = np.load('cvs/cvs_{}.npy'.format(channel))
    dvs = np.load('dvs/dvs_{}.npy'.format(channel))
    avs = np.load('avs/avs_{}.npy'.format(channel))
    counts = np.load('counts/counts_{}.npy'.format(channel))

    er_c = np.load('er_c/er_c_{}.npy'.format(channel))
    er_d = np.load('er_d/er_d_{}.npy'.format(channel))
    er_a = np.load('er_a/er_a_{}.npy'.format(channel))

    plt.errorbar(cvs,avs,xerr=er_c,yerr=er_a,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("c vs a: Dataset 4(19/6/22),  $-cP_1(t)+(d)(P_2(t)+2)$",fontdict={'fontsize': 18})
plt.xlabel('c',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('a',fontsize=18)
plt.savefig("c_a.png")



 
plt.figure(figsize = (16,9))
for channel in range(4):
    
    cvs = np.load('cvs/cvs_{}.npy'.format(channel))
    dvs = np.load('dvs/dvs_{}.npy'.format(channel))
    avs = np.load('avs/avs_{}.npy'.format(channel))
    counts = np.load('counts/counts_{}.npy'.format(channel))

    er_c = np.load('er_c/er_c_{}.npy'.format(channel))
    er_d = np.load('er_d/er_d_{}.npy'.format(channel))
    er_a = np.load('er_a/er_a_{}.npy'.format(channel))

    plt.errorbar(avs,dvs,xerr=er_a,yerr=er_d,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("a vs d: Dataset 4(19/6/22),  $-cP_1(t)+(d)(P_2(t)+2)$",fontdict={'fontsize': 18})
plt.xlabel('a',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('d',fontsize=18)
plt.savefig("a_d.png")
