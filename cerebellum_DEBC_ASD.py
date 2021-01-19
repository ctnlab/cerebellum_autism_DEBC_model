'''
@Author Piernadrea Mirino,  Daniele Caligiore, Emiliano Trimarco  2020-2021
'''
import nest as nt
import random
import numpy as np
import math as m
import datetime
from os.path import exists
import os


delay_F = round(random.uniform(50., 50.), 2)
delay_E = round(random.uniform(50., 50.), 2)

# Kernel function slow
def kernel_slow (x) :
    

    a=15
    b=15
    c=0.75
    d=0.2
    e=1.8
    f=1.3
    
    try:
        y=a*m.exp(-(m.fabs((x+d)*b)**e)/f)*(-m.sin((x+d)/m.exp(1)))**c  
    except:   
        y=0      

    return y

# Kernel function fast
def kernel_fast (x) :
    
    a=15
    b=15
    c=0.75
    d=0.1
    e=1.8
    f=1.3
    
    try:
        y=a*m.exp(-(m.fabs((x+d)*b)**e)/f)*(-m.sin((x+d)/m.exp(1)))**c
    except:
        y=0

    return y


#calculate time within range -400, 0 (IO time)
def window(TIO, t, wnd=0.40):

    t1= (t - (TIO - wnd))-wnd

    return t1


# spike generator function, pass arrays with values (time) that you want to use
def SpikeGenerator(sequence_spike_time=[100.]):

    spike = nt.Create("spike_generator", 1, {"spike_times": sequence_spike_time, 'precise_times': True})

    return spike


# Write  file
def scriviFile (obj, path=os.getcwd()+"/mice/", output_file="outputKernel.txt"):
    
    if not os.path.isdir(path):
         os.mkdir(path)
    file = path+output_file
                
    if (exists(file)):
        out_file = open(file, "a")
        out_file.write("%s\n" % str(obj))
        out_file.close()
    else:
        out_file = open(file, "w")
        out_file.write("%s\n" % str(obj))
        out_file.close()

#firing rate
def firingRate(spike_train,sig=100):
     kernel_gau=el.kernels.GaussianKernel(sigma=sig*pq.ms)

     firing_rate=el.statistics.instantaneous_rate(neo.SpikeTrain(spike_train, units=pq.ms, t_start=np.array(spike_train).min(), t_stop=np.array(spike_train).max()), sampling_period=pq.Quantity([1], 'ms'), kernel=kernel_gau)
    
     return firing_rate
       
                    
# This function run the two ways cerebellum network and return an complex object
def runNetwork(wf=0.5, ws=0.5, wpfc=0.5, seed_mice=0, w_inc_pfc_i=0.5, weight_gain_i=30, w_io_bf_i=500., w_io_bs_i=500., w_cs_stim=0.5, w_inc_pfc_i_nsp=0.00001, weight_k_increment_i_s=0.00001, weight_k_increment_i_f=0.00001, Ie_pfc_random = 370., Ie_m1_random = 370.,  w_n_nd_f_to_n_m1 = 100., w_n_nd_s_to_n_pfc=100.): 
    #==============================================================================
    #     ISI - interstimulus interval - the time distance in ms between GC spikes and IO spikes
    #     wf - weight fast line
    #     ws- weight slow line
    #     wpfc- weight pfc to m1
    #     seed_mice - seed
    #     w_inc_pfc_i - weigth incrementi pfc to m1
    #     w_inc_pfc_i_nsp=  weigth incrementi pfc to m1 no spike
    #     weight_k_increment_i_s weigth increment gc -> pj slow
    #     weight_k_increment_i_f weigth increment gc -> pj fast
    #     weight_gain_i- correction kernel
    #     w_io_bf_i weight io on purkinje fast line
    #     w_io_bs_i weight io on purkinje slow line 
    #     w_cs_stim weight Mossy fiber(stim) to Dentate Nucleus
    #==============================================================================
  
    #==============================================================================
    # Setup Parameters 
    #==============================================================================
    
    
    print ("EXP ---- W_FAST ----- >>>>> " + str(wf))
    print ("EXP ---- W_SLOW ----- >>>>> " + str(ws))
    print ("EXP ---- W_PFC ----- >>>>> " + str(wpfc))
        
    
    nt.ResetKernel() # reset nest
    nt.SetKernelStatus({'grng_seed': seed_mice})
    nt.SetKernelStatus({"local_num_threads": 4})
                
    #number of neuron test                                                     
    nNeuronGR_Fast = 768
    nNeuronGR_Slow = 768
    
    nNeuronPJ_Fast = 15
    nNeuronPJ_Slow = 15   

    nNeuronND_Fast = 12
    nNeuronND_Slow = 12

    nNeuronIO_Fast = 15
    nNeuronIO_Slow = 15

    nNeuronM1 = 12
    nNeuronPFC = 12
    
               
    # neuron type
    neuron_Type= 'iaf_psc_alpha'

    ISI = -20
    wnd=0.40 # kernel window
    version = 10    
    start_network = datetime.datetime.today() #start elaboration time
   
    #loop parameter
    dcfrom = 0 #start
    dcstep = 10 #rate
    dcto = 2500 #end
    
    weight_k_increment_pfc_to_m1 = w_inc_pfc_i #costant increment weight pfc to m1
    weight_k_increment_pfc_to_m1_no_spike =  w_inc_pfc_i_nsp
    weight_k_increment_s = weight_k_increment_i_s #costant increment weight Granular cells to PJ
    weight_k_increment_f = weight_k_increment_i_f #costant increment weight Granular cells to PJ
    
    weight_gain= weight_gain_i #normalization parallel fiber
    sum_w_fast=sum_w_slow=sum_w_pfc= 0
    
    # Initialization w
    wi_af_bf = wf
    wi_as_bs = ws
    wi_pfc_m1 = wpfc                
    w_io_bf = w_io_bf_i
    w_io_bs = w_io_bs_i
    
    #variables
    w_fast=[]
    w_slow=[]
    t_chart_fast= []
    w_chart_fast= []
    t_chart_slow= []
    w_chart_slow= []
    t_chart_pfc= []
    w_chart_pfc= []
    w_pfc_to_m1=[]
    w_pfc_to_m1_no_spike=[]
    buffer_kernel_k_fast=[]
    buffer_kernel_t_fast=[]
    buffer_kernel_k_slow=[]
    buffer_kernel_t_slow=[]
    spike_pfc_count=[]
    spike_io_count_f = 0
    spike_io_count_s = 0
    spike_pj_count_f=0
    spike_pj_count_s=0
    spike_pfc_count=0
    spike_m1_count=0    
     
    def generateSpikeTrainIO(step=10, start=1300+ISI, dStim=20, isiStim=780, sTrain=2500):         
         # Train of impulse on IO  
        s_io= []
        for i in range(sTrain/isiStim) :    
             s_io += SpikeGenerator(map(float,(range(start, start+dStim, step))))
             start += isiStim+dStim    
        return s_io
         
    
    def generateSpikeTrainGR(step=25, start=1000, dStim=300, isiStim=500, sTrain=2500):         
         # Train of impulse on Granular Cells CS         
        s_gr= []
        for i in range(sTrain/isiStim) :    
             s_gr += SpikeGenerator(map(float,(range(start, start+dStim, step))))
             start += isiStim+dStim    
        return s_gr
        
    s_io = generateSpikeTrainIO()
    
    s_gc = generateSpikeTrainGR()
    
    #single shot CS
    s_gc_base = generateSpikeTrainGR(25, 100, 300, 500, 600)
     
    noise_in = nt.Create("poisson_generator")
    nt.SetStatus(noise_in, {"rate": 2500.0, "start": 1., "stop": 2500.})
    
   
    # Spike detectors
    spike_det_gc_f = nt.Create("spike_detector")
    spike_det_gc_s = nt.Create("spike_detector")
    spike_det_nb_f = nt.Create("spike_detector")
    spike_det_nb_s = nt.Create("spike_detector")    
    spike_det_nio_f = nt.Create("spike_detector")
    spike_det_nio_s = nt.Create("spike_detector")
    spike_det_m1 = nt.Create("spike_detector")
    spike_det_pfc = nt.Create("spike_detector")
    spike_det_ndf = nt.Create("spike_detector")####
    spike_det_nds = nt.Create("spike_detector")####
    
    # Voltmeters granular cells
    vm_gc_f = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_gc_f, {'n_events': 0})
    vm_gc_s = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_gc_s, {'n_events': 0})

    # Voltmeters purkinje cells
    vm_pjf = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_pjf, {'n_events': 0})

    vm_pjL = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_pjL, {'n_events': 0})
    
    #==============================================================================
    # Create nodes 
    #==============================================================================
    
    # Granular Cells
    n_a_f = nt.Create(neuron_Type, nNeuronGR_Fast)
    n_a_s = nt.Create(neuron_Type, nNeuronGR_Slow)
    nt.SetStatus(n_a_f,  {'I_e': 370.0, "tau_m": 10., "V_th": -55.0})
    nt.SetStatus(n_a_s,  {'I_e': 370.0, "tau_m": 10., "V_th": -55.0})

    #Purkinje
    n_b_f = nt.Create(neuron_Type, nNeuronPJ_Fast)
    n_b_s = nt.Create(neuron_Type, nNeuronPJ_Slow)
    nt.SetStatus(n_b_f,  {'I_e': 380.0, "tau_m": 10.0, "V_th": -55.0})	
    nt.SetStatus(n_b_s,  {'I_e': 380.0, "tau_m": 10.0, "V_th": -55.0})
    
    # Inferior Olive
    n_io_f = nt.Create(neuron_Type, nNeuronIO_Fast)
    n_io_s = nt.Create(neuron_Type, nNeuronIO_Slow)
    nt.SetStatus(n_io_f,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    nt.SetStatus(n_io_s,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    
    # Dentate Nucleus
    n_nd_f = nt.Create(neuron_Type, nNeuronND_Fast)
    n_nd_s = nt.Create(neuron_Type, nNeuronND_Slow)
    nt.SetStatus(n_nd_f,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    nt.SetStatus(n_nd_s,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    
    # Motor primary cortex
    n_m1 = nt.Create(neuron_Type, nNeuronM1)
    nt.SetStatus(n_m1,  {'I_e': Ie_m1_random, "tau_m": 10. , "V_th": -55.0})
    
    # prefrontal cortex
    n_pfc = nt.Create(neuron_Type, nNeuronPFC)
    nt.SetStatus(n_pfc,  {'I_e': Ie_pfc_random , "tau_m": 10., "V_th": -55.0})
    
    
    #==============================================================================
    # Connections
    #==============================================================================
      
    nt.Connect(s_gc, n_a_f, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": delay_E}) #CS to GCf
    nt.Connect(s_gc, n_a_s, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": delay_E}) #CS to GCs
    nt.Connect(s_gc, n_nd_f, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": delay_E})  #CS to ND
    nt.Connect(s_gc, n_nd_s, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": delay_E})  #CS to ND

    nt.Connect(s_gc_base, n_a_f, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": delay_E}) #CS to GCf
    nt.Connect(s_gc_base, n_a_s, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": delay_E}) #CS to GCs   
    nt.Connect(s_gc_base, n_nd_f, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": delay_E})  #CS to ND
    nt.Connect(s_gc_base, n_nd_s, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": delay_E})  #CS to ND
    
    nt.Connect(s_io, n_io_f, syn_spec={'model': 'static_synapse', 'weight': 100., "delay": delay_E}) #US to IO
    nt.Connect(s_io, n_io_s, syn_spec={'model': 'static_synapse', 'weight': 100., "delay": delay_E}) #US to IO
    nt.Connect(n_io_f, n_b_f, 'one_to_one', syn_spec={'model': 'static_synapse', 'weight': w_io_bf}) #IO to PJ
    nt.Connect(n_io_s, n_b_s, 'one_to_one', syn_spec={'model': 'static_synapse', 'weight': w_io_bs}) #IO to PJ 
    nt.Connect(n_io_f, n_nd_f, syn_spec={'model': 'static_synapse', 'weight': 60.}) #Io to ND 
    nt.Connect(n_io_s, n_nd_s, syn_spec={'model': 'static_synapse', 'weight': 60.}) #Io to ND 
    
    # Connect PJ to ND divide by 2
    nt.Connect(n_b_f, n_nd_f, syn_spec={'model': 'static_synapse', 'weight': -7})
    nt.Connect(n_b_s, n_nd_s, syn_spec={'model': 'static_synapse', 'weight': -7})  
    
    nt.Connect(n_nd_f, n_m1, syn_spec={'model': 'static_synapse', 'weight': w_n_nd_f_to_n_m1, "delay": delay_F }) #w dentate nucleus to m1 and pfc
    nt.Connect(n_nd_s, n_pfc, syn_spec={'model': 'static_synapse', 'weight': w_n_nd_s_to_n_pfc}) #w dentate nucleus to m1 and pfc
    nt.Connect(noise_in, n_nd_f, syn_spec={'model': 'static_synapse', 'weight': round(random.uniform(0.1, 0.5), 2)}) #w dentate nucleus to m1 and pfc	
    nt.Connect(noise_in, n_nd_s, syn_spec={'model': 'static_synapse', 'weight': round(random.uniform(0.1, 0.5), 2)}) #w dentate nucleus to m1 and pfc
 
    
    # Connect voltimeters to neurons
    nt.Connect(vm_gc_f, n_a_f, 'all_to_all')
    nt.Connect(vm_gc_s, n_a_s, 'all_to_all')
    nt.Connect(vm_pjf, n_b_f, 'all_to_all') #purkinie line fast
    nt.Connect(vm_pjL, n_b_s, 'all_to_all') #purkinie line slow

    # Connect spike detector to neurons   
    nt.Connect(n_a_f, spike_det_gc_f)    
    nt.Connect(n_a_s, spike_det_gc_s)
    nt.Connect(n_b_f, spike_det_nb_f)
    nt.Connect(n_b_s, spike_det_nb_s)
    nt.Connect(n_io_f, spike_det_nio_f)
    nt.Connect(n_io_s, spike_det_nio_s)
    nt.Connect(n_m1, spike_det_m1)
    nt.Connect(n_pfc, spike_det_pfc)
    nt.Connect(n_nd_f, spike_det_ndf)######
    nt.Connect(n_nd_s, spike_det_nds)######
    
    #==============================================================================
    # Loop
    #==============================================================================
    for i, sim in enumerate(range(dcfrom, dcto, dcstep)):
    
        # Membrane potential and time Granular Cells Fast
        t_gc_f = nt.GetStatus(vm_gc_f, 'events')[0]['times']
        v_gc_f = nt.GetStatus(vm_gc_s, 'events')[0]['V_m']
        
        # Membrane potential and time Granular Cells Slow
        t_gc_s = nt.GetStatus(vm_gc_f, 'events')[0]['times']
        v_gc_s = nt.GetStatus(vm_gc_s, 'events')[0]['V_m']
        
        # Membrane potential and time putkinje line fast
        t_pjf = nt.GetStatus(vm_pjf, 'events')[0]['times']
        v_pjf = nt.GetStatus(vm_pjf, 'events')[0]['V_m']
        
        # Membrane potential and time putkinje line fast
        t_pjL = nt.GetStatus(vm_pjL, 'events')[0]['times']
        v_pjL = nt.GetStatus(vm_pjL, 'events')[0]['V_m']
                    
        dSD_gc_f = nt.GetStatus(spike_det_gc_f, keys="events")[0]
        evs_gc_f = dSD_gc_f["senders"]
        ts_gc_f = dSD_gc_f["times"]
 
        dSD_gc_s = nt.GetStatus(spike_det_gc_s, keys="events")[0]
        evs_gc_s = dSD_gc_s["senders"]
        ts_gc_s = dSD_gc_s["times"]       
            
        # Spike Purkinje Cells
        dSD_b_f = nt.GetStatus(spike_det_nb_f, keys="events")[0]
        evs_b_f = dSD_b_f["senders"]
        ts_b_f = dSD_b_f["times"]    
                                
        dSD_b_s = nt.GetStatus(spike_det_nb_s, keys="events")[0]
        evs_b_s = dSD_b_s["senders"]
        ts_b_s = dSD_b_s["times"]
    
        # Spike Inferior Olive
        dSD_io_f = nt.GetStatus(spike_det_nio_f, keys="events")[0]
        evs_io_f = dSD_io_f["senders"]
        ts_io_f = dSD_io_f["times"]
    
        # Spike Inferior Olive
        dSD_io_s = nt.GetStatus(spike_det_nio_s, keys="events")[0]
        evs_io_s = dSD_io_s["senders"]
        ts_io_s = dSD_io_s["times"]
  
    
        # Spike M1
        dSD_m1 = nt.GetStatus(spike_det_m1, keys="events")[0]
        evs_m1 = dSD_m1["senders"]
        ts_m1 = dSD_m1["times"]
    
        # Spike pfC
        dSD_pfc = nt.GetStatus(spike_det_pfc, keys="events")[0]
        evs_pfc = dSD_pfc["senders"]
        ts_pfc = dSD_pfc["times"]

        # Spike ndf
        dSD_ndf = nt.GetStatus(spike_det_ndf, keys="events")[0]
        evs_ndf = dSD_ndf["senders"]
        ts_ndf = dSD_ndf["times"]

        # Spike nds
        dSD_nds = nt.GetStatus(spike_det_nds, keys="events")[0]
        evs_nds = dSD_nds["senders"]
        ts_nds = dSD_nds["times"]
        

        ####################################################
        #################### LINE FAST ####################
        ####################################################
    
        if len(evs_gc_f) >=1 :
            #print ("LF1")

     
            # Record all time Spike of Granular Cells
            tm_f = np.unique(ts_gc_f*0.001)
        
            # Time Spike Inferior Olive e Purkinje
            TIO_f =  (float(nt.GetStatus(n_io_f, 't_spike')[0])*0.001)
            TPJ_f =  (float(nt.GetStatus(n_b_f, 't_spike')[0])*0.001) 
                                
                                    
            # If inferior olive send a Spike
            if  np.count_nonzero(evs_io_f) > spike_io_count_f:
                
                # Current number spike of inferior olive 
                spike_io_count_f = np.count_nonzero(evs_io_f)
                
    
                # For each Granular spike in the interval between TIO and TIO-tm_i
                for tm_i in tm_f:
             
                    if (TIO_f - tm_i) <= wnd :
                       
                        # Calculate the time in the range -800 / 0
                        window_tm_f = window(TIO_f, tm_i)

                        # Create an array of weights calculated with the kernel fast
                        w_fast.append(kernel_fast(window_tm_f))
                                  
                        # Parameters for the chart fast
                        buffer_kernel_k_fast.append(kernel_fast(window_tm_f))
                        buffer_kernel_t_fast.append(window_tm_f)
    
    
            for tm_i in tm_f:
            
               # If pj f send a Spike
               if  np.count_nonzero(evs_b_f) > spike_pj_count_f:
    
                   # Current number spike of pj f
                   spike_pj_count_f = np.count_nonzero(evs_b_f)
     
                   if (m.fabs(TPJ_f-tm_i) > wnd and (sum_w_fast <= wi_af_bf)) : w_fast.append(-weight_k_increment_f)                          
          
                 



        ####################################################
        #################### LINE SLOW ####################
        ####################################################
    
        if len(evs_gc_s) >=1 :
                                      
            # Record all time Spike of Granular Cells
            tm_s = np.unique(ts_gc_s*0.001)
        
            # Time Spike Inferior Olive e Purkinje
            TIO_s =  (float(nt.GetStatus(n_io_s, 't_spike')[0])*0.001)
            TPJ_s =  (float(nt.GetStatus(n_b_s, 't_spike')[0])*0.001) 
                                    
                                    
            # If inferior olive send a Spike
            if  np.count_nonzero(evs_io_s) > spike_io_count_s:
                
                # Current number spike of inferior olive 
                spike_io_count_s = np.count_nonzero(evs_io_s)
                
    
                # For each Granular spike in the interval between TIO and TIO-tm_i
                for tm_i in tm_s:
                    
        
                    if (TIO_s - tm_i) <= wnd :
                        
        
                        # Calculate the time in the range
                        window_tm_s = window(TIO_s, tm_i)
                                                             
                        # Create an array of weights calculated with the kernel slow
                        w_slow.append(kernel_slow(window_tm_s))
                                    
                        # Parameters for the chart slow
                        buffer_kernel_k_slow.append(kernel_slow(window_tm_s))
                        buffer_kernel_t_slow.append(window_tm_s)
    
            for tm_i in tm_s:
    
               # If pj s send a Spike
               if  np.count_nonzero(evs_b_s) > spike_pj_count_s:
    
                   # Current number spike of pj s
                   spike_pj_count_s = np.count_nonzero(evs_b_s)
                
                   if (m.fabs(TPJ_s-tm_i) > wnd and (sum_w_slow <= wi_as_bs)) : w_slow.append(-weight_k_increment_s)                               
                    

              
        #if sum_w_fast <= wi_af_bf  : w_fast.append(-weight_k_increment) # Costant inccreaase weight if sum_w_fast <= wi_af_bf
        sum_w_fast = wi_af_bf - (sum(w_fast) *  weight_gain)

        
    
        #if sum_w_slow <= wi_as_bs : w_slow.append(-weight_k_increment) #costant increaase weight
        sum_w_slow = wi_as_bs - (sum(w_slow) * weight_gain)
        
        
                                
        #weight chart
        t_chart_fast.append(i)
        w_chart_fast.append(sum_w_fast)
        t_chart_slow.append(i)
        w_chart_slow.append(sum_w_slow)
                                
        # If pfc  send a Spike
        TM1 = np.array(nt.GetStatus(spike_det_m1, keys="events")[0]["times"])[-1:]*0.001   
        TPFC = np.array(nt.GetStatus(spike_det_pfc, keys="events")[0]["times"])[-1:]*0.001      
    
        if (np.count_nonzero(evs_pfc) > spike_pfc_count) and (np.count_nonzero(evs_m1) > spike_m1_count) and m.fabs(TM1-TPFC) < 0.40:
             spike_pfc_count = np.count_nonzero(evs_pfc)
             spike_m1_count = np.count_nonzero(evs_m1)                     
             w_pfc_to_m1.append(weight_k_increment_pfc_to_m1) #incrementare weight on array pfc 
        else:
            if  np.count_nonzero(evs_m1) > spike_m1_count:
                w_pfc_to_m1_no_spike.append(weight_k_increment_pfc_to_m1_no_spike)
                spike_m1_count = np.count_nonzero(evs_m1)
            
            
        sum_w_pfc=wi_pfc_m1 + sum(w_pfc_to_m1) - sum(w_pfc_to_m1_no_spike) # Sum array weigth
                                
        t_chart_pfc.append(i)
        w_chart_pfc.append(sum_w_pfc)
    
        syn_dict_fast= {"model": "static_synapse", "weight": sum_w_fast}
        syn_dict_slow = {"model": "static_synapse", "weight": sum_w_slow}
        syn_dict_pfc_to_m1 = {"model": "static_synapse", "weight": sum_w_pfc} #peso tra pfc e m1
    
        #  Connect GRf to PJs
        nt.Connect(n_a_f, n_b_f, syn_spec=syn_dict_fast)
        # Connect GRs to PJs
        nt.Connect(n_a_s, n_b_s, syn_spec=syn_dict_slow)
        # Connect pfc to m1
        nt.Connect(n_pfc, n_m1, syn_spec=syn_dict_pfc_to_m1)
            
        nt.Simulate(dcstep) #12/01 inserito la variabile step invece del valore fisso 1
    
     
    
    obj_returned={"timesSpikeGC_f":ts_gc_f,"sendersGC_f":evs_gc_f,
                  "timesSpikeGC_s":ts_gc_s,"sendersGC_s":evs_gc_s,
                  "potentialGC_f":np.array(v_gc_f, dtype=float),
                  "tPotentialGC_f":np.array(t_gc_f, dtype=float),
                  "potentialGC_s":np.array(v_gc_s, dtype=float),
                  "tPotentialGC_s":np.array(t_gc_s, dtype=float),                  
                  "timesSpikePjF":ts_b_f,"sendersPjF":evs_b_f, 
                  "timesSpikePjS": ts_b_s, "sendersPjS":evs_b_s, 
                  "timesSpikeIO_f": ts_io_f, "sendersIO_f":evs_io_f,
                  "timesSpikeIO_s": ts_io_s, "sendersIO_s":evs_io_s, 
                  "timesSpikeM1": ts_m1, "timesSpikePFC": ts_pfc,
                  "timesSpikeNDf": ts_ndf, "timesSpikeNDs": ts_nds, 
                  "timesSpikeStartNetwork": start_network,
                  "tw_chart_fast": t_chart_fast, "w_chart_fast": w_chart_fast,
                  "tw_chart_slow": t_chart_slow, "w_chart_slow": w_chart_slow,
                  "t_chart_pfc": t_chart_pfc, "w_chart_pfc": w_chart_pfc, 
                  "senderM1": evs_m1, "senderPFC": evs_pfc,
                  "potentialPjA":np.array(v_pjf, dtype=float),
                  "tPotentialPjA":np.array(t_pjf, dtype=float),
                  "potentialPjL":np.array(v_pjL, dtype=float),
                  "tPotentialPjL":np.array(t_pjL, dtype=float)}
                  

    return  obj_returned

      
#test = runNetwork ()
