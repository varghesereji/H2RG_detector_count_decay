import numpy as np
import matplotlib.pyplot as plt

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']
plt.figure(figsize = (16,9))
for channel in range(4):
    
    pvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/pvs/pvs_{}.npy'.format(channel))
    qvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/qvs/qvs_{}.npy'.format(channel))
    counts = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/counts/counts_{}.npy'.format(channel))

    er_p = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/er_p/er_p_{}.npy'.format(channel))
    er_q = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/er_q/er_q_{}.npy'.format(channel))

    plt.errorbar(counts,pvs,yerr=er_p,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("p vs count: Dataset 4(19/6/22),  $px+q(-\dfrac{3}{8}x+\dfrac{1}{2}x^2)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('p',fontsize=18)
plt.savefig("p_legendre.png")




plt.figure(figsize = (16,9))
for channel in range(4):
    
    
    pvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/pvs/pvs_{}.npy'.format(channel))
    qvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/qvs/qvs_{}.npy'.format(channel))
    counts = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/counts/counts_{}.npy'.format(channel))

    er_p = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/er_p/er_p_{}.npy'.format(channel))
    er_q = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/er_q/er_q_{}.npy'.format(channel))

    plt.errorbar(counts,qvs,yerr=er_q,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("q vs count: Dataset 4(19/6/22),  $px+q(-\dfrac{3}{8}x+\dfrac{1}{2}x^2)$",fontdict={'fontsize': 18})
plt.xlabel('Count(ADU)',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('q',fontsize=18)
plt.savefig("q_legendre.png")



plt.figure(figsize = (16,9))
for channel in range(4):
    
    
    pvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/pvs/pvs_{}.npy'.format(channel))
    qvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/qvs/qvs_{}.npy'.format(channel))
    counts = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/counts/counts_{}.npy'.format(channel))

    er_p = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/er_p/er_p_{}.npy'.format(channel))
    er_q = np.load('/home/varghese/Desktop/DS5_DP1/Codes/9_Orthogonalize/New_orthogonal/er_q/er_q_{}.npy'.format(channel))

    plt.errorbar(pvs,qvs,xerr=er_p,yerr=er_q,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))


plt.title("p vs q: Dataset 4(19/6/22),  $px+q(-\dfrac{3}{8}x+\dfrac{1}{2}x^2)$",fontdict={'fontsize': 18})
plt.xlabel('p',fontsize=18)
plt.legend()
#plt.ylim(0.0,0.65)
plt.ylabel('q',fontsize=18)
plt.savefig("p_q.png")
