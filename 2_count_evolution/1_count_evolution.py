import numpy as np
import matplotlib.pyplot as plt



Data_0G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G/AvgDataCube_0G\
_2.npy')

Data_1G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_1G/AvgDataCube_1G\
_2.npy')

Data_2G = np.load('/home/varghese/Desktop/DS5_DP1/Dataset/Avg_2G/AvgDataCube_2G\
_2.npy')


mask = (Data_0G[-1] > 20000) & (Data_0G[-1] < 30000)


plt.figure(figsize=(16,9))



plt.plot(Data_2G[:,mask],color='blue',alpha=0.4)
plt.plot(Data_1G[:,mask],color='green',alpha=0.4)
plt.plot(Data_0G[:,mask],color='red',alpha=0.4)

plt.savefig('Count_evolution_DS5.png')
