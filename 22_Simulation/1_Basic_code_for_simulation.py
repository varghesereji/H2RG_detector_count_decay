import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

L, M, N = 1.8354467816716038, 0.23986204637532926, 0.0980528855003108

def exponential_from_p(t,p):
    q = L*p**2+M*p+N
    c = p+(3/8)*q
    b=q/c
    a= c**2/q
    return a*(np.exp(-b*t/20)-1)

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


c_0 = 0
dt = 1  # time is in the unit of per readout
rate = 3000 * dt  # Rate is in the unit of per readout.
t_max = 10 / dt

up_the_ramp = generate_up_the_ramp(c_0,  rate,  t_max,  dt)

pvs = np.load('/home/varghese/Desktop/DS5_DP1/Codes/22_Simulation/Parameter_values/new_pvs_0.npy')
count = np.load('/home/varghese/Desktop/DS5_DP1/Codes/22_Simulation/Parameter_values/new_count_0.npy')

interp_p = interp1d(count, pvs, fill_value='extrapolate')

updated_count = up_the_ramp[-1]
updated_p = interp_p(updated_count)
time_up_the_ramp = np.linspace(0,t_max,len(up_the_ramp))
time_decay = np.arange(t_max, 30/dt ,dt)
decay_curve = updated_count \
    + (updated_count / 100 ) \
    * ( exponential_from_p(time_decay-time_decay[0], updated_p))

plt.figure()
plt.plot(time_up_the_ramp, up_the_ramp,'o-')
plt.plot(time_decay ,decay_curve,'o-')
plt.show(block=False)

print('Difference is:',33000 - up_the_ramp[-1])
