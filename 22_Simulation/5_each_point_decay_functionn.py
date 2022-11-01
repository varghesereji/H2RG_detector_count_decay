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
    L, M, N = 1.8354467816716038, 0.23986204637532926, 0.098052885500310
    q = L*p**2+M*p+N
    c = p+(3/8)*q
    b=q/c
    a= c**2/q
    return (a/100) * (np.exp(-b*t/20)-1)

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

def generate_up_the_ramp_decay_added(c_0, rate, t_max, dt):
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
        updated_p = interpolating_0(c_0, channel)
        updated_count = c_t + exponential_from_p(dt,updated_p) * c_t
        print(c_t,exponential_from_p(dt,updated_p) * c_t,updated_count)
        up_the_ramp_values.append(c_t)
        c_0 = updated_count
        t = t + dt
    return up_the_ramp_values

def generate_simulation_each_instant_decay(c_0=0, rate=3000,t_max=10, dt=1, channel=0):
    rate = rate * dt
    t_max = 10 / dt
    
    up_the_ramp = generate_up_the_ramp_decay_added(c_0,  rate,  t_max,  dt, channel)
    updated_count = up_the_ramp[-1]
    updated_p = interp_p(updated_count)
    time_up_the_ramp = np.linspace(0,t_max,len(up_the_ramp))
    time_decay = np.arange(t_max, 30/dt ,dt)
    decay_curve = updated_count \
        + (updated_count) \
        * ( exponential_from_p(time_decay-time_decay[0], updated_p))
    simulated_curve = np.hstack((up_the_ramp, decay_curve))
    simulated_time = np.hstack((time_up_the_ramp.time_decay))
    return simulated_curve, simulated_time 
