import sys
import HxRGproc.reduction
from glob import glob
from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt



path= '/home/varghese/Desktop/20220619/hpf'

order_list = sorted(os.listdir(path))



order_list_firsthalf = order_list[2:]


AvgDataCube = None

l=0

for directory in order_list_firsthalf:
       
    file_path=os.path.join(path,directory,'fits')
    #print(file_path)
    print('Called files in', file_path)
    file_list=sorted(os.listdir(file_path))
    #print(file_path,file_list)

    if file_list==[]:
        pass
    else:
        list_with_path = []
        for files in file_list:
            UTRlistT = os.path.join(file_path,files)

            list_with_path.append(UTRlistT)

        UTRfiles = sorted(list_with_path, key=HxRGproc.reduction.instruments.sort_filename_key_function_HPFLinux)

        RawDataCube = np.array([fits.getdata(utr) for utr in UTRfiles])
        DataCube = HxRGproc.reduction.remove_biases_in_cube(RawDataCube,no_channels=4, do_LSQmedian_correction=-99999)

        
        if l%3 ==0:
            np.save('/home/varghese/Desktop/DS5_DP1/Dataset/1G_dataset_all/Exposure_{}.npy'.format(l//3),DataCube)
            print('/home/varghese/Desktop/DS5_DP1/Dataset/1G_dataset_all/Exposure_{}.npy Saved'.format(l//3))
        elif l%3 == 1:
            np.save('/home/varghese/Desktop/DS5_DP1/Dataset/2G_dataset_all/Exposure_{}.npy'.format(l//3),DataCube)
            print('/home/varghese/Desktop/DS5_DP1/Dataset/2G_dataset_all/Exposure_{}.npy Saved'.format(l//3))
        else:
            np.save('/home/varghese/Desktop/DS5_DP1/Dataset/0G_dataset_all/Exposure_{}.npy'.format(l//3),DataCube)
            print('/home/varghese/Desktop/DS5_DP1/Dataset/0G_dataset_all/Exposure_{}.npy Saved'.format(l//3))
      
    l+=1
