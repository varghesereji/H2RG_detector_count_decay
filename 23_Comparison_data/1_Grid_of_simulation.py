import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn

from decay_simulations import *

c_rate = np.arange(1000,4500,1)
difference_grid = None
t_max = 26
dt = 1
c_rate_axis = np.arange(1000,4500,500)
for rate in c_rate:
    print(rate)
    ideal_ramp = generate_up_the_ramp(0,rate,30,dt)
    simulated_ramp, simulation_time = generate_simulation_each_instant_decay(0, rate, t_max, dt, 0)
    #print(ideal_ramp, simulated_ramp)
    Relative_difference = (ideal_ramp - simulated_ramp) * 100 / ideal_ramp
    #print(Relative_difference[1:20])
    if difference_grid is None:
        difference_grid = Relative_difference[1:t_max-1]
    else:
        difference_grid = np.vstack((difference_grid,Relative_difference[1:t_max-1]))
    #for num in range(len(difference_grid)):
        #print(ideal_ramp[num],simulated_ramp[num],Relative_difference[num])
    # fig, axs = plt.subplots(2)
    # fig.suptitle('rate = {}'.format(rate))
    # axs[0].plot(ideal_ramp[1:t_max-1],'o-')
    # axs[0].plot(simulated_ramp[1:t_max-1],'o-')
    # axs[1].plot(Relative_difference[1:t_max-1],'o-')
    # plt.show(block=False)
#print(difference_grid)
plt.figure()
#plt.imshow(difference_grid)
mapping = sn.heatmap(data=difference_grid, xticklabels=np.arange(1,t_max-1,1), yticklabels=c_rate)
plt.title('Percentage change')
plt.savefig('Decay_amount_s1.png')
#plt.show(block=False)
