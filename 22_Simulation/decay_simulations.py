import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def interpolating_p(updated_count,channel):
    pvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/22_Simulation/Parameter_values/new_pvs_{}.npy'.format(channel))
    count = np.load('/home/varghese/Desktop/DS5_DP1/Codes/22_Simulation/Parameter_values/new_count_{}.npy'.format(channel))
    interp_p = interp1d(count, pvs, fill_value='extrapolate')
    updated_p = interp_p(updated_count)
    return updated_p

def exponential_from_p(t,p):
    L, M, N = 1.8354467816716038, 0.23986204637532926, 0.0980528855003108
    q = L*p**2+M*p+N
    c = p+(3/8)*q
    b=q/c
    a= c**2/q
    return (a/100)*(np.exp(-b*t/20)-1)

def find_c_t(c_0, rate, dt):
    '''
    c_0 : the initial value
    rate : Rate of increase in count per unit time
    dt : infititesimal time.
    '''
    c_t = c_0 + rate*dt
    return c_t

def generate_up_the_ramp(c_0, rate, t_max, dt):
    '''
    c_0 : Initial value
    t_max : Maximum time to continue this process
    rate : Rate of increase in count per unit time
    dt : infititesimal time.
    '''
    up_the_ramp_values = [c_0]
    t = 0
    while t<=t_max:
        c_t = find_c_t(c_0,  rate, dt)
        up_the_ramp_values.append(c_t)
        c_0=c_t
        t=t+dt
    return up_the_ramp_values

def generate_decay_curve(count,time,dt,channel):
    '''
    Parameters
    -----------
    count : value of count at that instant
    time : decay time of the particular instant
    dt : infinitesimal time

    Return
    --------------
    decay_curve: array of functional values of decay
    '''
    time_decay = np.arange(time, 30/dt ,dt)
    #print('count_1',count)
    updated_p = interpolating_p(count, channel)
    ##print(updated_p)
    decay_curve = (count) \
                * ( exponential_from_p(time_decay-time_decay[0], updated_p))
    return decay_curve

'''
------------------------------------------------------------------------
-----------Expected overestimation--------------------------------------
------------------------------------------------------------------------
'''

def generate_up_the_ramp_decay_added(c_0, rate, t_max, dt, channel):
    '''
    c_0 : Initial value
    t_max : Maximum time to continue this process
    rate : Rate of increase in count per unit time
    dt : infititesimal time.

    Return
    ---------------------------------------
    up_the_ramp_values: list. Each element is the sum of increased count and decay because of previous instant.
    '''
    up_the_ramp_values = [c_0]
    t = 0
    while t<=t_max:
        c_t = find_c_t(c_0,  rate, dt)
        updated_p = interpolating_p(c_0, channel)
        updated_count = c_t + exponential_from_p(dt,updated_p) * c_t
       # print(c_t,exponential_from_p(dt,updated_p) * c_t,updated_count)
        up_the_ramp_values.append(c_t)
        c_0 = updated_count
        t = t + dt
    return up_the_ramp_values

def generate_simulation_each_instant_decay(c_0=0, rate=3000,t_max=10, dt=1, channel=0):
    rate = rate * dt
    t_max = 10 / dt
    
    up_the_ramp = generate_up_the_ramp_decay_added(c_0,  rate,  t_max,  dt, channel)
    updated_count = up_the_ramp[-1]
    updated_p = interpolating_p(updated_count, channel)
    time_up_the_ramp = np.linspace(0,t_max,len(up_the_ramp))
    time_decay = np.arange(t_max, 30/dt ,dt)
    decay_curve = updated_count \
        + (updated_count) \
        * ( exponential_from_p(time_decay-time_decay[0], updated_p)) #Decau after closing of shutter
    simulated_curve = np.hstack((up_the_ramp, decay_curve))
    simulated_time = np.hstack((time_up_the_ramp, time_decay))
    return simulated_curve, simulated_time 


'''
------------------------------------------------------------------------
-----------Expected underestimation--------------------------------------
------------------------------------------------------------------------
'''


def adding_old_new_decay(old, new):
    decay_sum = old[1:] + new
    return decay_sum

def decay_added_up_the_ramp(c_0, rate, t_max, dt, channel):
    up_the_ramp_values = [c_0]
    t = 0
    while t<=t_max:
        c_t = find_c_t(c_0,rate,dt)
        ##print('c_t',c_t)
        exp_decay = generate_decay_curve(c_t, t, dt, channel)
        if t==0:
            old_decay = exp_decay
        else:
            old_decay = adding_old_new_decay(old_decay,exp_decay)
        updated_c = c_t + old_decay[0]
        ##print(c_t,old_decay[0],updated_c)
        up_the_ramp_values.append(updated_c)
        t +=dt
        c_0 = updated_c
    return up_the_ramp_values

def generate_simulation_previous_decay_sum(c_0=0, rate=3000, t_max=10, dt=1, channel=0):
    '''
    Parameters
    ---------------------------------
    c_0 : initial value of count
    rate : increase in count after each readout (Default value=3000)
    t_max : time which shutter is closing. In the unit of per readout ( Default value = 10)
    dt : time difference between two successive readout (Default value = 1)
    channel : Number of channel. (Default value = 0)
    Returns
    ---------------------------------
    simulated_curve : the simulated curve which includs sum of value of decay formed from previous instants
    simulation_time : time of simulation 
    '''
    rate = rate * dt
    t_max = 10 / dt
    ##print('c_o', c_0)
    up_the_ramp = decay_added_up_the_ramp(c_0,  rate,  t_max,  dt, channel)
    time_up_the_ramp = np.linspace(0,t_max,len(up_the_ramp))
    decay_curve = generate_decay_curve(up_the_ramp[-1],t_max,dt, channel) + up_the_ramp[-1]
    time_decay = np.arange(t_max,30/dt,dt)
    simulated_curve = np.hstack((up_the_ramp, decay_curve))
    simulation_time = np.hstack((time_up_the_ramp, time_decay))
    return simulated_curve, simulation_time
