

import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import os
import multiprocessing
import subprocess
from matplotlib.backends.backend_pdf import PdfPages

color_list = ['black','blue','green','red','dimgray','deepskyblue','lightgreen','lightcoral']

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
    V1 = v1(t)
    V2 = v2(t)
    return -p*V1+q*V2

def Constants(l,m,n):  # In terms of p
    def exp_p(t,p):
        #q = (l*p**2+m*p+n)
        #c = (p+(3/8)*(l*p**2+m*p+n))
        #b= ((l*p**2+m*p+n) / (p+(3/8)*(l*p**2+m*p+n)))
        #a= ((p+(3/8)*(l*p**2+m*p+n))**2/ (l*p**2+m*p+n))
        return ((p+(3/8)*(l*p**2+m*p+n))**2/ (l*p**2+m*p+n))*(np.exp(-((l*p**2+m*p+n) / (p+(3/8)*(l*p**2+m*p+n)))*t)-1)
    return exp_p

def Constants1(l,m,n): # In terms of q
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

def Fitting_1(path,name):
    pvs = []
    qvs = []
    er_p = []
    er_q = []
    counts = []
    for channel in range(4):
        Filename = name+'{}.npy'.format(channel)
        File = os.path.join(path,Filename)
        dc = np.load(File)
        Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
        Sigma1 = np.median(np.sqrt(Sigma2.flatten()))+np.zeros(20)
        fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:512*(channel+1)]
        interval = 250
        for flux_bin in range(500,60000,interval):
            pix_mask = (dc[-1] > (flux_bin-interval//2)) & (dc[-1] < (flux_bin+interval//2)) & (fiber_mask == 1)
            avg_count = np.nanmedian(dc[10,pix_mask])
            delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
            time1 = np.arange(0,np.size(delta_change_cube,axis=0),1) / (np.size(delta_change_cube,axis=0))
            if not np.isnan(avg_count) and delta_change_cube[-1]<0:
                Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
                try:
                    fitting_data_variables, error = curve_fit(polynomial,time1,delta_change_cube,np.array([1.0,0.5]),sigma = Sigma)
                    p,q = fitting_data_variables
                    dp = np.sqrt(error[0][0])
                    dq = np.sqrt(error[1][1])
                    pvs.append(p)
                    qvs.append(q)
                    counts.append(flux_bin)
                    er_p.append(dp)
                    er_q.append(dq)
                except Exception as e:
                    print(e)
                    print('{0},{1}'.format(avg_count,channel),'cannot be done')
    print("channel", channel, 'done pvs and qvs')
    return pvs,qvs,counts,er_p,er_q




