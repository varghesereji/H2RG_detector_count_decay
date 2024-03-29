

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time
import multiprocessing

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

t1 = time.time()

def objective(x,l,m,n):
    return l*x**2+m*x+n

def exponential(t,a,b):
    return a*(np.exp(-b*t)-1)

def ab_from_pq(p,q):
    c = p+(3/8)*q
    d=q

    b=d/c
    a=c/b
    return a,b


def v1(x):
    return x

def v2(x):
    return ((-3/8)*x+x**2/2)


def polynomial(t,p,q):
#    print(t)                                                                                                                                                                                              
    V1 = v1(t)
    V2 = v2(t)

    return -p*V1+q*V2

def Constants(l,m,n):
    
    def exp_p(t,p):
        #q = (l*p**2+m*p+n)
        #c = (p+(3/8)*(l*p**2+m*p+n))
        #b= ((l*p**2+m*p+n) / (p+(3/8)*(l*p**2+m*p+n)))
        #a= ((p+(3/8)*(l*p**2+m*p+n))**2/ (l*p**2+m*p+n))

        return ((p+(3/8)*(l*p**2+m*p+n))**2/ (l*p**2+m*p+n))*(np.exp(-((l*p**2+m*p+n) / (p+(3/8)*(l*p**2+m*p+n)))*t)-1)
    return exp_p

def Constants1(l,m,n):
    
    def exp_p(t,p):
        # c = (((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q)
        # d = q
        # b = d/c = (q/(((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))
        # a = c**2 / d = ((((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))**2/(q)
        return ((((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))**2/(q)*(np.exp(-(q/(((-m+np.sqrt(m**2-4*l*(n-q)))/(2*l))+(3/8)*q))*t)-1)
    return exp_p


def curve_from_p(t,p,l,m,n):
    q = l*p**2+m*p+n
    c = p+(3/8)*q
    b=q/c
    a= c**2/q

    return a*(np.exp(-b*t)-1)


    
def Fitting(channel):
    pvs=np.load('pvs/pvs_{}.npy'.format(channel))
    qvs=np.load('qvs/qvs_{}.npy'.format(channel))
    er_p=np.load('er_p/er_p_{}.npy'.format(channel))
    er_q=np.load('er_q/er_q_{}.npy'.format(channel))
    avs=np.load('avs/avs_{}.npy'.format(channel))
    bvs=np.load('bvs/bvs_{}.npy'.format(channel))
    er_a=np.load('er_a/er_a_{}.npy'.format(channel))
    er_b=np.load('er_b/er_b_{}.npy'.format(channel))

    counts = np.load('counts/counts_{}.npy'.format(channel))
    
    a_n, b_n = ab_from_pq(pvs,qvs)

    fitted_values,error = curve_fit(objective,pvs,qvs,np.array([0.5,0.5,0.0]),sigma=er_q)
    
    l,m,w = fitted_values
    er_l = np.sqrt(error[0][0])
    er_m = np.sqrt(error[1][1])
    er_e=np.sqrt(error[2][2])
    #print(channel,l,m,e,error)
    a_n, b_n = ab_from_pq(pvs,objective(pvs,l,m,w))
    # Here, p is q actually
    p_values = []
    p_err = []
    Count_new = []
    dc = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G_{}.npy'.format(channel))
    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)
    interval = 100
    flag = 0
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    for flux_bin in range(500,60000,interval):
#rows = []
        
        
        pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc[10,pix_mask])
        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1)
        time1=time1/time1[-1]
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            P_ini=0.05
            flag+=1
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))

            try:
                PQ,er_pq = curve_fit(polynomial,time1,delta_change_cube,np.array([0.05,0.05]),sigma=Sigma)
                P_new,Q_new = PQ

                ab, er_ab = curve_fit(exponential,time1,delta_change_cube,np.array([0.1,0.01]),sigma=Sigma)
                a,b = ab
                                                                                 
                P_values, error = curve_fit(Constants(l,m,w),time1,delta_change_cube,np.array([P_ini]),sigma=Sigma)
                P = P_values[0]
                perr =  np.sqrt(error[0][0])
                print('channel',channel,'p found',P,'initial p',P_ini,'error',perr)
                p_values.append(P)
                p_err.append(perr)
                Count_new.append(flux_bin)
                plt.figure(figsize=(16,9))
                plt.plot(time1,delta_change_cube,'o-',label='data')
                plt.plot(time1,exponential(time1,a,b),'--',label='exponential from a and b')   # Polynomial with a and b
                plt.plot(time1,curve_from_p(time1,P,l,m,w),'--',label='exp finding p={}'.format(P))
                #plt.plot(time1,curve_from_p(time1,P_ini,l,m,w),'--',label='exp model from initial value of p={}'.format(P_ini))
                plt.plot(time1,polynomial(time1,P_new,Q_new),'--',label='polynomial with p={0},q={0}'.format(P_new,Q_new))   # Polynomial with p and q
                plt.plot(time1,polynomial(time1,P_new,objective(P_new,l,m,w)),'--',label='polynomial with p only'.format(P_new)) # Polynomial with p only
                plt.xlabel('time')
                plt.legend()
                plt.ylabel('count')
                plt.title('p Fitting Channel:{0}, Count{1}'.format(channel,flux_bin))
                plt.savefig('p_model/p_model_{0}_{1}.png'.format(channel,flux_bin))
            except Exception as e:
                print(e)
    print(channel,'done')
    np.save('new_count_{}.npy'.format(channel),Count_new)
    np.save('new_p_err_{}.npy'.format(channel),p_err)
    np.save('new_pvs_{}.npy'.format(channel),p_values)
    
'''
    for n in range(len(counts)):
        c = counts[n]
        delta_change_cube=np.load('Ch_{0}/Decay_curve_{1}_{2}.npy'.format(channel,channel,c))
        Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
        Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)
        Sigma = Sigma1*100/(c*np.sqrt(1000))
        
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / (np.size(delta_change_cube,axis=0))
        P_values, error = curve_fit(Constants(l,m,n),time1,delta_change_cube,np.array([0.4]),sigma=Sigma)

        print(P_values,error)
        plt.figure(figsize=(16,9))
        plt.plot(time1,delta_change_cube,label='Data')
        plt.plot(time1,exponential(time1,avs[n],bvs[n]),label='Exponential_fit')
        plt.plot(time1,polynomial(time1,pvs[n],qvs[n]),label='Polynomian fit')
        plt.plot(time1,polynomial(time1,pvs[n],objective(pvs[n],l,m,e)),label='Polynomial from model')
        
        plt.plot(time1,exponential(time1,a_n[n],b_n[n]),label='Generated from model')
        
        plt.legend()
        plt.title('count={6},a={0},b={1},p={2},q={3},a_n={4},b_n={5}'.format(avs[n],bvs[n],pvs[n],qvs[n],a_n[n],b_n[n],c))
        plt.savefig('Ch_{0}/Model_comparison_{1}_{2}.png'.format(channel,channel,c))

        p_values.append(P_values[0])
    
    #plt.figure()
    #plt.errorbar(pvs,qvs,xerr=er_p,yerr=er_q,fmt='.',color=color_list[channel],ecolor=color_list[4+channel],elinewidth=1,capsize=10,label='Channel {}'.format(channel+1))
    #plt.plot(pvs,objective(pvs,l,m,n),'--',color=color_list[channel])
    #plt.savefig("Fitted_model_{}.png".format(channel))
    #print("Saved")

    np.save('pvalues_{}.npy'.format(channel),p_values)
    print(channel,'done')
    plt.figure()
    plt.plot(counts,p_values)
    plt.savefig('p_{}.png'.format(channel))
''' 

process1 = multiprocessing.Process(target=Fitting,args = [0])
process2 = multiprocessing.Process(target=Fitting,args = [1])
process3 = multiprocessing.Process(target=Fitting,args = [2])
process4 = multiprocessing.Process(target=Fitting,args = [3])


process1.start()
process2.start()
process3.start()
process4.start()

process1.join()
process2.join()
process3.join()
process4.join()


