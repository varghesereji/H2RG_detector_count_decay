import matplotlib.pyplot as plt

from decay_simulations import *
from modeling_all_available_datacubes import *





def Fitting_2(channel,path,name):
    
    pvs,qvs,counts,er_p,er_q = Fitting_1(path,name)    
    fitted_values,error = curve_fit(objective,pvs,qvs,np.array([0.5,0.5,0.0]),sigma=er_q)

    l,m,w = fitted_values
    er_l = np.sqrt(error[0][0])
    er_m = np.sqrt(error[1][1])
    er_e=np.sqrt(error[2][2])
    print(channel,'l:',l,'m:',m,'w:',w,'errors',er_l,er_m,er_e)                                                                                                                                      
    #a_n, b_n = ab_from_pq(pvs,objective(pvs,l,m,w))
    # Here, p is q actually                                                                                                                                          
    p_values = []
    p_err = []
    Count_new = []

    Filename = name+'{}.npy'.format(channel)
    File = os.path.join(path,Filename)
    dc = np.load(File)
    pdfplots = PdfPages("Decay_curves_ch_{0}.pdf".format(channel))
    Sigma2 = np.load('/home/varghese/Desktop/DP1/Codes_Dataset_3/9_readout_error/Long_exp_Variance_channel_{}.npy'.format(channel))
    Sigma1 = np.median(np.sqrt(Sigma2.flatten())) + np.zeros(20)
    interval = 250
    flag = 0
    fiber_mask = np.load('/home/varghese/Desktop/DP1/Masks/Needed_fiber_mask.npy')[:,512*channel:(512*(channel+1))]
    for flux_bin in range(500,60000,interval):
        pix_mask = (dc[-1] > (flux_bin-interval/2)) & (dc[-1] < (flux_bin+interval/2)) & (fiber_mask ==1)
        avg_count = np.nanmedian(dc[10,pix_mask])
        delta_change_cube = np.nanmedian((dc[10:,pix_mask]-dc[10,pix_mask])*100/dc[10,pix_mask],axis=1)
        time1 = np.arange(0,np.size(delta_change_cube,axis=0),1)
        time1=time1/time1[-1]
        rate = avg_count / 11
        c_0 = np.nanmedian(dc[0])
        t_max = 9
        dt = 1
        each_instant_decay_added, simulation_time_2 = generate_simulation_each_instant_decay(c_0, rate, t_max, dt, channel)
        dc_instant_decay = (each_instant_decay_added[t_max+3:] - each_instant_decay_added[t_max+3])*100/each_instant_decay_added[t_max+3]
        previous_decay_sum,simulation_time_1 = generate_simulation_previous_decay_sum(c_0, rate, t_max, dt, channel)
        dc_decay_sum = (previous_decay_sum[t_max+3:] - previous_decay_sum[t_max+3])*100/previous_decay_sum[t_max+3]
        if not np.isnan(avg_count) and delta_change_cube[-1]<0:
            P_ini=0.05
            flag+=1
            Sigma = Sigma1*100/(avg_count*np.sqrt(np.sum(pix_mask)))
            try:
                PQ,er_pq = curve_fit(polynomial,time1,delta_change_cube,np.array([0.05,0.05]),sigma=Sigma)
                P_new,Q_new = PQ
                ab, er_ab = curve_fit(exponential,time1,delta_change_cube,np.array([0.1,0.01]),sigma=Sigma)
                a,b = ab
                P_values, error = curve_fit(Constants(l,m,w),time1,delta_change_cube,np.array([P_ini]),sigma=Sigma)
                P = P_values[0]
                perr =  np.sqrt(error[0][0])
                p_values.append(P)
                p_err.append(perr)
                Count_new.append(flux_bin)
                
                fig = plt.figure(figsize=(16,9))
                plt.plot(time1,delta_change_cube,'o-',label='data')
                plt.plot(time1,polynomial(time1,P_new,objective(P_new,l,m,w)),'--',label='p polynomial'.format(P_new))
                plt.plot(time1, dc_instant_decay,'o-',label='instant decay')
                plt.plot(time1, dc_decay_sum,'o-',label='sum of decays')
                plt.xlabel('time')
                plt.legend()
                plt.ylabel('count')
                plt.title('p Fitting Channel:{0}, Count{1}'.format(channel,flux_bin))
                pdfplots.savefig()
            except Exception as e:
                print(e)
    print(channel,'done pvs')
    pdfplots.close()
    #PVs[channel] = p_values
    #Counts[channel] = Count_new
    #er_P[channel] = p_err

#Plotting

path = '/home/varghese/Desktop/DS5_DP1/Dataset/Avg_0G'
File = 'AvgDataCube_0G_'

#DS5: AvgDataCube_0G_

process1 = multiprocessing.Process(target=Fitting_2,args = [0,path,File])
process2 = multiprocessing.Process(target=Fitting_2,args = [1,path,File])
process3 = multiprocessing.Process(target=Fitting_2,args = [2,path,File])
process4 = multiprocessing.Process(target=Fitting_2,args = [3,path,File])

process1.start()
process2.start()
process3.start()
process4.start()

process1.join()
process2.join()
process3.join()
process4.join()
