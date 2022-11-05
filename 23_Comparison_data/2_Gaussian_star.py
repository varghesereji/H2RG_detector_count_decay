import numpy as np
import matplotlib.pyplot as plt

from decay_simulations import *

def fake_star(x):
    return 30000 * np.exp(-x**2)

X = np.arange(-2,2,0.1)

Ideal_star = fake_star(X)

t_max = 20

count_rate = Ideal_star / t_max

simulated_star = []
for rate in count_rate:
    simulated_ramp, simulation_time = generate_simulation_each_instant_decay(0,rate,t_max,1,0)
    simulated_star.append(simulated_ramp[t_max-1])

simulated_star = np.array(simulated_star)
plt.figure()
plt.plot(X,Ideal_star, label='Ideal star')
plt.plot(X,simulated_star, label='simulated_star')
plt.legend()
plt.savefig('simulated_star.png')
