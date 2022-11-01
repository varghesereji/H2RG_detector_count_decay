import numpy as np


DataCube = ['/home/varghese/Desktop/DS5_DP1/Dataset/0G_dataset_all/Exposure_{}.npy'.format(file) for file in range(16,20,1)]

AvgDataCube = None

l = 0
for order,file in enumerate(DataCube):
    dc = np.load(file)

    if AvgDataCube is None:
        AvgDataCube = dc
    else:
        AvgDataCube = AvgDataCube + dc

    print(l,'done')
    l=l+1

AvgDataCube = AvgDataCube / l

AvgFluxPerFrame = np.mean(np.diff(AvgDataCube[1:9],axis=0))
AvgDataCube = AvgDataCube+AvgFluxPerFrame

for i,subcube in enumerate(np.array_split(AvgDataCube,4,axis=2)):                                                                                     
    np.save('/home/varghese/Desktop/DS5_DP1/Codes/19_Checking_trapping/AvgDataCube_l4_0G_{}.npy'.format(i),subcube)
    print('/home/varghese/Desktop/DS5_DP1/Codes/19_Checking_trapping/AvgDataCube_l4_0G_{}.npy Saved'.format(i))

