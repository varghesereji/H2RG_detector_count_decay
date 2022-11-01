from decay_simulations import *
import matplotlib.pyplot as plt


previous_decay_sum,simulation_time_1 = generate_simulation_previous_decay_sum()

each_instant_decay_added, simulation_time_2 = generate_simulation_each_instant_decay()

plt.figure()
plt.plot(simulation_time_1,previous_decay_sum,'o-',label='prev decay sum')
plt.plot(simulation_time_2,each_instant_decay_added,'o-',label='instant decay')
plt.ylim(33000-700,33000)
plt.legend()
plt.savefig('compare_simulation.png')
