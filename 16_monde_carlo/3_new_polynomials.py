

import numpy as np
from scipy.optimize import curve_fit
import multiprocessing
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from scipy.optimize import least_squares

t1 = time.time()

def v1(x):
    return x+1
def v2(x):
    return (3/4)*(2*x**2+x-1)

def objective(t,b,c):
    V1=v1(t)
    V2=v2(t)
    return b*V1+c*V2 


channel = 1

dc = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))
Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
Sigma1 = np.median(np.sqrt(Sigma2.flatten()))
#fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512(channel+1))]

interval = 100
flux_bin = 20000

pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2))# & (fiber_mask == 1)

avg_count = np.nanmedian(dc[10,pix_mask])
delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
Sigma_array = np.zeros(20)+Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

Sigma1 = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
time = []
Fake_data = None
N = 1000
for t in range(0,20):
    time.append(t)
    mean = delta_change_cube[t]
    Sigma = Sigma1
    if t == 0:
        Fake_data = np.zeros(N)
    else:
        D_point = np.random.normal(mean,Sigma,N)
        Fake_data=np.vstack((Fake_data,D_point))

Fake_data = np.transpose(Fake_data)
time1 = np.array(time)/10-1

avs = []
bvs = []
cvs = []
er_a = []
er_b = []
er_c = []
flag = 0
for data in Fake_data:
    print(data)
    try:
        print(flag)
    
        fitting_data_variables, error = curve_fit(objective,time1,data,np.array([0.5,0.5]),sigma=Sigma_array)
        print(time1)
        b,c = fitting_data_variables
        #da = np.sqrt(error[0][0])
        db = np.sqrt(error[0][0])
        dc1 = np.sqrt(error[1][1])
        print(error)
        plt.figure()
        plt.plot(time1,data,'o-')
        plt.plot(time1,objective(time1,b,c),"--")
        plt.plot(time1,delta_change_cube,">-")
        plt.savefig('Plot_{}_new.png'.format(flag))
        flag+=1
        #avs.append(a)
        bvs.append(b)
        cvs.append(c)
        
        #er_a.append(da)
        er_b.append(db)
        er_c.append(dc)
    except Exception as e:
        print(e)
'''
plt.figure()

print(avs)
#plt.errorbar(avs,bvs,xerr=er_a,yerr=er_b,fmt='.',elinewidth=1,capsize=10)
plt.plot(avs,bvs,'.')
plt.xlabel("a")
plt.ylabel("b")
plt.savefig("a_vs_b.png")
''




plt.figure()

print(avs)
#plt.errorbar(avs,bvs,xerr=er_a,yerr=er_b,fmt='.',elinewidth=1,capsize=10)
plt.plot(avs,cvs,'.')
plt.xlabel("a")
plt.ylabel("c")
plt.savefig("a_vs_c.png")

'''
plt.figure()

#print(avs)
#plt.errorbar(avs,bvs,xerr=er_a,yerr=er_b,fmt='.',elinewidth=1,capsize=10)
plt.plot(cvs,bvs,'.')
plt.xlabel("c")
plt.ylabel("b")
plt.savefig("c_vs_b_new_polynomials.png")
