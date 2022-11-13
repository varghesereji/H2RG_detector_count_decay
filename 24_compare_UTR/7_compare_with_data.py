import matplotlib.pyplot as plt

from decay_simulations import *
from modeling_all_available_datacubes import *


def Fitting_2(channel,path,name):
    
    w = 0
    Filename = name+'{}.npy'.format(channel)
    File = os.path.join(path,Filename)
    dc = np.load(File)
    pdfplots = PdfPages("Decay_curves_ch_{0}.pdf".format(channel))
    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)
    interval = 100
    flag = 0
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    for flux_bin in range(500,60000,interval):
        print(channel, w)
        w = w + 1
        pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc[10,pix_mask])
        UTR_data = np.mean(dc[:,pix_mask],axis=1)
        c_0 = np.nanmedian(dc[0])
        t_max = 8
        dt = 1
        rate = np.mean(UTR_data[8]) / 9
        each_instant_decay_added, simulation_time_2 = generate_simulation_each_instant_decay(c_0, rate, t_max, dt, channel)
        previous_decay_sum,simulation_time_1 = generate_simulation_previous_decay_sum(c_0, rate, t_max, dt, channel)
        if not np.isnan(avg_count):
            try:
                fig = plt.figure(figsize=(16,9))
                plt.plot(UTR_data[:9], 'o-', label='data')
                plt.plot(each_instant_decay_added[:9],'o-',label='instant decay')
                plt.plot(previous_decay_sum[:9],'o-',label='sum of decays')
                plt.xlabel('time')
                plt.legend()
                plt.ylabel('count')
                plt.title('p Fitting Channel:{0}, Count{1}'.format(channel,flux_bin))
                pdfplots.savefig()
            except Exception as e:
                print(e)
    print(channel,'done pvs')
    pdfplots.close()

path = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G'
File = 'AvgDataCube_0G_'

# process1 = multiprocessing.Process(target=Fitting_2,args = [0,path,File])
# process2 = multiprocessing.Process(target=Fitting_2,args = [1,path,File])
# process3 = multiprocessing.Process(target=Fitting_2,args = [2,path,File])
# process4 = multiprocessing.Process(target=Fitting_2,args = [3,path,File])

# process1.start()
# process2.start()
# process3.start()
# process4.start()

# process1.join()
# process2.join()
# process3.join()
# process4.join()

for q in range(4):
    Fitting_2(q, path, File)
