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


# Kernel function cognitive pathway
def kernel_cognitive (x) :
    

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

# Kernel function motor pathway
def kernel_motor (x) :
    
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
def runNetwork(wm=0.5, wc=0.5, wpfc=0.5, seed_mice=0, w_inc_pfc_i=0.5, weight_gain_i=30, w_io_bm_i=500., w_io_bc_i=500., w_cs_stim=0.5, w_inc_pfc_i_nsp=0.00001, weight_k_increment_i_c=0.00001, weight_k_increment_i_m=0.00001, Ie_pfc_random = 370., Ie_m1_random = 370.,  w_n_nd_m_to_n_m1 = 100., w_n_nd_c_to_n_pfc=100.): 
    #==============================================================================
    #     ISI - interstimulus interval - the time distance in ms between GC spikes and IO spikes
    #     wm - weight motor line
    #     wc- weight cognitive line
    #     wpfc- weight pfc to m1
    #     seed_mice - seed
    #     w_inc_pfc_i - weigth incrementi pfc to m1
    #     w_inc_pfc_i_nsp=  weigth incrementi pfc to m1 no spike
    #     weight_k_increment_i_c weigth increment gc -> pj cognitve
    #     weight_k_increment_i_m weigth increment gc -> pj motor
    #     weight_gain_i- correction kernel
    #     w_io_bm_i weight io on purkinje motor line
    #     w_io_bc_i weight io on purkinje cognitive line 
    #     w_cs_stim weight Mossy fiber(stim) to Dentate Nucleus
    #==============================================================================
  
    #==============================================================================
    # Setup Parameters 
    #==============================================================================
    
    
    print ("EXP ---- W_MOTOR ----- >>>>> " + str(wm))
    print ("EXP ---- W_COGNITIVE ----- >>>>> " + str(wc))
    print ("EXP ---- W_PFC ----- >>>>> " + str(wpfc))
        
    
    nt.ResetKernel() # reset nest
    nt.SetKernelStatus({'grng_seed': seed_mice})
    nt.SetKernelStatus({"local_num_threads": 4})
                
    #number of neuron test                                                     
    nNeuronGR_Motor = 768
    nNeuronGR_Cognitive = 768
    
    nNeuronPJ_Motor = 24
    nNeuronPJ_Cognitive = 24   

    nNeuronND_Motor = 12
    nNeuronND_Cognitive = 12

    nNeuronIO_Motor = 24
    nNeuronIO_Cognitive = 24

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
    weight_k_increment_c = weight_k_increment_i_c #costant increment weight Granular cells to PJ
    weight_k_increment_m = weight_k_increment_i_m #costant increment weight Granular cells to PJ
    
    weight_gain= weight_gain_i #normalization parallel fiber
    sum_w_motor=sum_w_cognitive=sum_w_pfc= 0
    
    # Initialization w
    wi_am_bm = wm
    wi_ac_bc = wc
    wi_pfc_m1 = wpfc                
    w_io_bm = w_io_bm_i
    w_io_bc = w_io_bc_i
    
    #variables
    w_motor=[]
    w_cognitive=[]
    t_chart_motor= []
    w_chart_motor= []
    t_chart_cognitive= []
    w_chart_cognitive= []
    t_chart_pfc= []
    w_chart_pfc= []
    w_pfc_to_m1=[]
    w_pfc_to_m1_no_spike=[]
    buffer_kernel_k_motor=[]
    buffer_kernel_t_motor=[]
    buffer_kernel_k_cognitive=[]
    buffer_kernel_t_cognitive=[]
    spike_pfc_count=[]
    spike_io_count_m = 0
    spike_io_count_c = 0
    spike_pj_count_m=0
    spike_pj_count_c=0
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
    spike_det_gc_m = nt.Create("spike_detector")
    spike_det_gc_c = nt.Create("spike_detector")
    spike_det_nb_m = nt.Create("spike_detector")
    spike_det_nb_c = nt.Create("spike_detector")    
    spike_det_nio_m = nt.Create("spike_detector")
    spike_det_nio_c = nt.Create("spike_detector")
    spike_det_m1 = nt.Create("spike_detector")
    spike_det_pfc = nt.Create("spike_detector")
    spike_det_nd_m = nt.Create("spike_detector")
    spike_det_nd_c = nt.Create("spike_detector")
    
    # Voltmeters granular cells
    vm_gc_m = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_gc_m, {'n_events': 0})
    vm_gc_c = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_gc_c, {'n_events': 0})

    # Voltmeters purkinje cells
    vm_pj_m = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_pj_m, {'n_events': 0})

    vm_pj_c = nt.Create('voltmeter', params={'to_file': False, 'withgid': True, 'withtime': True, 'interval': 0.1})
    nt.SetStatus(vm_pj_c, {'n_events': 0})
    
    #==============================================================================
    # Create nodes 
    #==============================================================================
    
    # Granular Cells
    n_a_m = nt.Create(neuron_Type, nNeuronGR_Motor)
    n_a_c = nt.Create(neuron_Type, nNeuronGR_Cognitive)
    nt.SetStatus(n_a_m,  {'I_e': 370.0, "tau_m": 10., "V_th": -55.0})
    nt.SetStatus(n_a_c,  {'I_e': 370.0, "tau_m": 10., "V_th": -55.0})

    #Purkinje
    n_b_m = nt.Create(neuron_Type, nNeuronPJ_Motor)
    n_b_c = nt.Create(neuron_Type, nNeuronPJ_Cognitive)
    nt.SetStatus(n_b_m,  {'I_e': 380.0, "tau_m": 10.0, "V_th": -55.0})	
    nt.SetStatus(n_b_c,  {'I_e': 380.0, "tau_m": 10.0, "V_th": -55.0})
    
    # Inferior Olive
    n_io_m = nt.Create(neuron_Type, nNeuronIO_Motor)
    n_io_c = nt.Create(neuron_Type, nNeuronIO_Cognitive)
    nt.SetStatus(n_io_m,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    nt.SetStatus(n_io_c,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    
    # Dentate Nucleus
    n_nd_m = nt.Create(neuron_Type, nNeuronND_Motor)
    n_nd_c = nt.Create(neuron_Type, nNeuronND_Cognitive)
    nt.SetStatus(n_nd_m,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    nt.SetStatus(n_nd_c,  {'I_e': 370.0, "tau_m": 10.0, "V_th": -55.0})
    
    # Motor primary cortex
    n_m1 = nt.Create(neuron_Type, nNeuronM1)
    nt.SetStatus(n_m1,  {'I_e': Ie_m1_random, "tau_m": 10. , "V_th": -55.0})
    
    # prefrontal cortex
    n_pfc = nt.Create(neuron_Type, nNeuronPFC)
    nt.SetStatus(n_pfc,  {'I_e': Ie_pfc_random , "tau_m": 10., "V_th": -55.0})
    
    
    #==============================================================================
    # Connections
    #==============================================================================
    
    #connection speed
    speed_stimuli = round(random.uniform(100., 100.), 2)# connection speed between sensory area and GR
    speed_dnm_m1 = round(random.uniform(100., 100.), 2)# connection speed between ND_m and M1
  
    nt.Connect(s_gc, n_a_m, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": speed_stimuli}) #CS to GC_m
    nt.Connect(s_gc, n_a_c, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": speed_stimuli}) #CS to GC_c
    nt.Connect(s_gc, n_nd_m, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": speed_stimuli})  #CS to ND
    nt.Connect(s_gc, n_nd_c, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": speed_stimuli})  #CS to ND

    nt.Connect(s_gc_base, n_a_m, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": speed_stimuli}) #CS to GC_m
    nt.Connect(s_gc_base, n_a_c, syn_spec={'model': 'static_synapse', 'weight': 500., "delay": speed_stimuli}) #CS to GC_c   
    nt.Connect(s_gc_base, n_nd_m, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": speed_stimuli})  #CS to ND
    nt.Connect(s_gc_base, n_nd_c, syn_spec={'model': 'static_synapse', 'weight': w_cs_stim, "delay": speed_stimuli})  #CS to ND
    
    nt.Connect(s_io, n_io_m, syn_spec={'model': 'static_synapse', 'weight': 100., "delay": speed_stimuli}) #US to IO
    nt.Connect(s_io, n_io_c, syn_spec={'model': 'static_synapse', 'weight': 100., "delay": speed_stimuli}) #US to IO
    nt.Connect(n_io_m, n_b_m, 'one_to_one', syn_spec={'model': 'static_synapse', 'weight': w_io_bm}) #IO to PJ
    nt.Connect(n_io_c, n_b_c, 'one_to_one', syn_spec={'model': 'static_synapse', 'weight': w_io_bc}) #IO to PJ 
    nt.Connect(n_io_m, n_nd_m, syn_spec={'model': 'static_synapse', 'weight': 60.}) #Io to ND 
    nt.Connect(n_io_c, n_nd_c, syn_spec={'model': 'static_synapse', 'weight': 60.}) #Io to ND 
    
    # Connect PJ to ND divide by 2
    nt.Connect(n_b_m, n_nd_m, syn_spec={'model': 'static_synapse', 'weight': -7})
    nt.Connect(n_b_c, n_nd_c, syn_spec={'model': 'static_synapse', 'weight': -7})  
    
    nt.Connect(n_nd_m, n_m1, syn_spec={'model': 'static_synapse', 'weight': w_n_nd_m_to_n_m1, "delay": speed_dnm_m1 }) #w dentate nucleus to m1 
    nt.Connect(n_nd_c, n_pfc, syn_spec={'model': 'static_synapse', 'weight': w_n_nd_c_to_n_pfc}) #w dentate nucleus to pfc
    nt.Connect(noise_in, n_nd_m, syn_spec={'model': 'static_synapse', 'weight': round(random.uniform(0.1, 0.5), 2)}) # dentate nucleus noise	
    nt.Connect(noise_in, n_nd_c, syn_spec={'model': 'static_synapse', 'weight': round(random.uniform(0.1, 0.5), 2)}) # dentate nucleus noise
 
    
    # Connect voltimeters to neurons
    nt.Connect(vm_gc_m, n_a_m, 'all_to_all')
    nt.Connect(vm_gc_c, n_a_c, 'all_to_all')
    nt.Connect(vm_pj_m, n_b_m, 'all_to_all') #purkinie line Motor
    nt.Connect(vm_pj_c, n_b_c, 'all_to_all') #purkinie line Cognitive

    # Connect spike detector to neurons   
    nt.Connect(n_a_m, spike_det_gc_m)    
    nt.Connect(n_a_c, spike_det_gc_c)
    nt.Connect(n_b_m, spike_det_nb_m)
    nt.Connect(n_b_c, spike_det_nb_c)
    nt.Connect(n_io_m, spike_det_nio_m)
    nt.Connect(n_io_c, spike_det_nio_c)
    nt.Connect(n_m1, spike_det_m1)
    nt.Connect(n_pfc, spike_det_pfc)
    nt.Connect(n_nd_m, spike_det_nd_m)
    nt.Connect(n_nd_c, spike_det_nd_c)
    
    #==============================================================================
    # Loop
    #==============================================================================
    for i, sim in enumerate(range(dcfrom, dcto, dcstep)):
    
        # Membrane potential and time Granular Cells Motor
        t_gc_m = nt.GetStatus(vm_gc_m, 'events')[0]['times']
        v_gc_m = nt.GetStatus(vm_gc_m, 'events')[0]['V_m']
        
        # Membrane potential and time Granular Cells Cognitive
        t_gc_c = nt.GetStatus(vm_gc_c, 'events')[0]['times']
        v_gc_c = nt.GetStatus(vm_gc_c, 'events')[0]['V_m']
        
        # Membrane potential and time putkinje line Motor
        t_pj_m = nt.GetStatus(vm_pj_m, 'events')[0]['times']
        v_pj_m = nt.GetStatus(vm_pj_m, 'events')[0]['V_m']
        
        # Membrane potential and time putkinje line Cognitive
        t_pj_c = nt.GetStatus(vm_pj_c, 'events')[0]['times']
        v_pj_c = nt.GetStatus(vm_pj_c, 'events')[0]['V_m']
                    
        dSD_gc_m = nt.GetStatus(spike_det_gc_m, keys="events")[0]
        evs_gc_m = dSD_gc_m["senders"]
        ts_gc_m = dSD_gc_m["times"]
 
        dSD_gc_c = nt.GetStatus(spike_det_gc_c, keys="events")[0]
        evs_gc_c = dSD_gc_c["senders"]
        ts_gc_c = dSD_gc_c["times"]       
            
        # Spike Purkinje Cells
        dSD_b_m = nt.GetStatus(spike_det_nb_m, keys="events")[0]
        evs_b_m = dSD_b_m["senders"]
        ts_b_m = dSD_b_m["times"]    
                                
        dSD_b_c = nt.GetStatus(spike_det_nb_c, keys="events")[0]
        evs_b_c = dSD_b_c["senders"]
        ts_b_c = dSD_b_c["times"]
    
        # Spike Inferior Olive
        dSD_io_m = nt.GetStatus(spike_det_nio_m, keys="events")[0]
        evs_io_m = dSD_io_m["senders"]
        ts_io_m = dSD_io_m["times"]
    
        # Spike Inferior Olive
        dSD_io_c = nt.GetStatus(spike_det_nio_c, keys="events")[0]
        evs_io_c = dSD_io_c["senders"]
        ts_io_c = dSD_io_c["times"]
  
    
        # Spike M1
        dSD_m1 = nt.GetStatus(spike_det_m1, keys="events")[0]
        evs_m1 = dSD_m1["senders"]
        ts_m1 = dSD_m1["times"]
    
        # Spike pfC
        dSD_pfc = nt.GetStatus(spike_det_pfc, keys="events")[0]
        evs_pfc = dSD_pfc["senders"]
        ts_pfc = dSD_pfc["times"]

        # Spike nd_m
        dSD_nd_m = nt.GetStatus(spike_det_nd_m, keys="events")[0]
        evs_nd_m = dSD_nd_m["senders"]
        ts_nd_m = dSD_nd_m["times"]

        # Spike nd_c
        dSD_nd_c = nt.GetStatus(spike_det_nd_c, keys="events")[0]
        evs_nd_c = dSD_nd_c["senders"]
        ts_nd_c = dSD_nd_c["times"]
        

        ####################################################
        #################### LINE MOTOR ####################
        ####################################################
    
        if len(evs_gc_m) >=1 :
            #print ("LF1")

     
            # Record all time Spike of Granular Cells
            tm_m = np.unique(ts_gc_m*0.001)
        
            # Time Spike Inferior Olive e Purkinje
            TIO_m =  (float(nt.GetStatus(n_io_m, 't_spike')[0])*0.001)
            TPJ_m =  (float(nt.GetStatus(n_b_m, 't_spike')[0])*0.001) 
                                
                                    
            # If inferior olive send a Spike
            if  np.count_nonzero(evs_io_m) > spike_io_count_m:
                
                # Current number spike of inferior olive 
                spike_io_count_m = np.count_nonzero(evs_io_m)
                
    
                # For each Granular spike in the interval between TIO and TIO-tm_i
                for tm_i in tm_m:
             
                    if (TIO_m - tm_i) <= wnd :
                       
                        # Calculate the time in the range -800 / 0
                        window_tm_m = window(TIO_m, tm_i)

                        # Create an array of weights calculated with the kernel motor
                        w_motor.append(kernel_motor(window_tm_m))
                                  
                        # Parameters for the chart motor
                        buffer_kernel_k_motor.append(kernel_motor(window_tm_m))
                        buffer_kernel_t_motor.append(window_tm_m)
    
    
            for tm_i in tm_m:
            
               # If pj_m send a Spike
               if  np.count_nonzero(evs_b_m) > spike_pj_count_m:
    
                   # Current number spike of pj_m
                   spike_pj_count_m = np.count_nonzero(evs_b_m)
     
                   if (m.fabs(TPJ_m-tm_i) > wnd and (sum_w_motor <= wi_am_bm)) : w_motor.append(-weight_k_increment_m)                          
          
                 



        ####################################################
        #################### LINE COGNITIVE ####################
        ####################################################
    
        if len(evs_gc_c) >=1 :
                                      
            # Record all time Spike of Granular Cells
            tm_c = np.unique(ts_gc_c*0.001)
        
            # Time Spike Inferior Olive e Purkinje
            TIO_c =  (float(nt.GetStatus(n_io_c, 't_spike')[0])*0.001)
            TPJ_c =  (float(nt.GetStatus(n_b_c, 't_spike')[0])*0.001) 
                                    
                                    
            # If inferior olive send a Spike
            if  np.count_nonzero(evs_io_c) > spike_io_count_c:
                
                # Current number spike of inferior olive 
                spike_io_count_c = np.count_nonzero(evs_io_c)
                
    
                # For each Granular spike in the interval between TIO and TIO-tm_i
                for tm_i in tm_c:
                    
        
                    if (TIO_c - tm_i) <= wnd :
                        
        
                        # Calculate the time in the range
                        window_tm_c = window(TIO_c, tm_i)
                                                             
                        # Create an array of weights calculated with the kernel cognitive
                        w_cognitive.append(kernel_cognitive(window_tm_c))
                                    
                        # Parameters for the chart cognitive
                        buffer_kernel_k_cognitive.append(kernel_cognitive(window_tm_c))
                        buffer_kernel_t_cognitive.append(window_tm_c)
    
            for tm_i in tm_c:
    
               # If pj_c send a Spike
               if  np.count_nonzero(evs_b_c) > spike_pj_count_c:
    
                   # Current number spike of pj_c
                   spike_pj_count_c = np.count_nonzero(evs_b_c)
                
                   if (m.fabs(TPJ_c-tm_i) > wnd and (sum_w_cognitive <= wi_ac_bc)) : w_cognitive.append(-weight_k_increment_c)                               
                    

              
        #if sum_w_motor <= wi_am_bm  : w_motor.append(-weight_k_increment) # Costant inccreaase weight if sum_w_motor <= wi_af_bf
        sum_w_motor = wi_am_bm - (sum(w_motor) *  weight_gain)

        
    
        #if sum_w_cognitive <= wi_ac_bc : w_cognitive.append(-weight_k_increment) #costant increaase weight
        sum_w_cognitive = wi_ac_bc - (sum(w_cognitive) * weight_gain)
        
        
                                
        #weight chart
        t_chart_motor.append(i)
        w_chart_motor.append(sum_w_motor)
        t_chart_cognitive.append(i)
        w_chart_cognitive.append(sum_w_cognitive)
                                
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
    
        syn_dict_motor= {"model": "static_synapse", "weight": sum_w_motor}
        syn_dict_cognitive = {"model": "static_synapse", "weight": sum_w_cognitive}
        syn_dict_pfc_to_m1 = {"model": "static_synapse", "weight": sum_w_pfc} #weight between pfc and m1
    
        #  Connect GR_m to PJ_m
        nt.Connect(n_a_m, n_b_m, syn_spec=syn_dict_motor)
        # Connect GR_c to PJ_c
        nt.Connect(n_a_c, n_b_c, syn_spec=syn_dict_cognitive)
        # Connect pfc to m1
        nt.Connect(n_pfc, n_m1, syn_spec=syn_dict_pfc_to_m1)
            
        nt.Simulate(dcstep) #
    
     
    
    obj_returned={"timesSpikeGC_m":ts_gc_m,"sendersGC_m":evs_gc_m,
                  "timesSpikeGC_c":ts_gc_c,"sendersGC_c":evs_gc_c,
                  "potentialGC_m":np.array(v_gc_m, dtype=float),
                  "tPotentialGC_m":np.array(t_gc_m, dtype=float),
                  "potentialGC_c":np.array(v_gc_c, dtype=float),
                  "tPotentialGC_c":np.array(t_gc_c, dtype=float),                  
                  "timesSpikePj_m":ts_b_m,"sendersPj_m":evs_b_m, 
                  "timesSpikePj_c": ts_b_c, "sendersPj_c":evs_b_c, 
                  "timesSpikeIO_m": ts_io_m, "sendersIO_m":evs_io_m,
                  "timesSpikeIO_c": ts_io_c, "sendersIO_c":evs_io_c, 
                  "timesSpikeM1": ts_m1, "timesSpikePFC": ts_pfc,
                  "timesSpikeND_m": ts_nd_m, "timesSpikeND_c": ts_nd_c, 
                  "timesSpikeStartNetwork": start_network,
                  "tw_chart_motor": t_chart_motor, "w_chart_motor": w_chart_motor,
                  "tw_chart_cognitive": t_chart_cognitive, "w_chart_cognitive": w_chart_cognitive,
                  "t_chart_pfc": t_chart_pfc, "w_chart_pfc": w_chart_pfc, 
                  "senderM1": evs_m1, "senderPFC": evs_pfc,
                  "potentialPj_m":np.array(v_pj_m, dtype=float),
                  "tPotentialPj_m":np.array(t_pj_m, dtype=float),
                  "potentialPj_c":np.array(v_pj_c, dtype=float),
                  "tPotentialPj_c":np.array(t_pj_c, dtype=float)}
                  

    return  obj_returned

      
#test = runNetwork ()
