'''
@Author Piernadrea Mirino,  Daniele Caligiore, Emiliano Trimarco  2020-2021
'''
import cerebellum_DEBC_ASD as mice
import numpy as np
import random
import os

#initial value of sections
number_of_mice = 15
number_of_section = 11

for n_mice in range(0, number_of_mice):
    
    
    seed_mice= random.randint(2500,3000)
    
    random.seed( seed_mice )
    
    #==============================================================================
    # set value of sections
    #==============================================================================
    
 
    wpfc = round(random.uniform(0.1, 0.1), 2) #w initializzation pfc     
    w_inc_pfc_i = round(random.uniform(0.2, 0.5), 2) #w inc pfc to m1 spike
    w_inc_pfc_i_nsp=round(random.uniform(0.015, 0.035), 3) #w pfc to m1 no spike 13/15
    
    wm = round(random.uniform(5., 5.), 2) #initializzation w kernel motor
    wc = round(random.uniform(20., 20.), 2) #initializzation w kernel cognitive 
    weight_k_increment_i_m=round(random.uniform(0.05, 0.05), 2) #0.05 costant w increment weight Granular cells to PJ no spike
    weight_k_increment_i_c=weight_k_increment_i_m

    w_cs_stim = round(random.uniform(500, 500), 2)#w CS stim 
    
    w_io_bm = random.randint(500, 500)*-1 #w US
    w_io_bm_i = w_io_bm #w io to pj  US
    w_io_bc_i = w_io_bm #w io to pj  US 
    
    weight_gain_i = round(random.uniform(5., 5.),2) #gain
    Ie_pfc_random = round(random.uniform(300.,365.),2)
    Ie_m1_random = round(random.uniform(300., 365.),2)
    
    w_n_nd_m_to_n_m1 = round(random.uniform(100., 100.), 2)
    w_n_nd_c_to_n_pfc = round(random.uniform(50., 50.), 2)
    
 
    path=os.getcwd()+"/session_mice_"+str(n_mice)+"/"
    
    mice.scriviFile(seed_mice, path, "seed.txt")
        
    ISI=20
  
    
    #==============================================================================
    # 
    #==============================================================================
    
    
    initial_value = "wm: " + str(wm) + "\n" + "wc: " + str(wc) + "\n" + "wpfc: " + str(wpfc) + "\n" + "seed_mice: "+ str(seed_mice) + "\n" + "w_inc_pfc_i: "+ str(w_inc_pfc_i) + "\n" + "w_inc_pfc_i_nsp: "+ str(w_inc_pfc_i_nsp) + "\n" + "weight_k_increment_i_m: "+ str(weight_k_increment_i_m) + "\n" + "weight_k_increment_i_c: "+ str(weight_k_increment_i_c) + "\n" + "weight_gain_i: " + str(+weight_gain_i) + "\n" + "w_io_bm_i: "+str(w_io_bm_i) + "\n" + "w_io_bc_i: "+ str(w_io_bc_i) + "\n" + "w_cs_stim: "+ str(w_cs_stim) + "\n" + "Ie_pfc_random: " + str(Ie_pfc_random) + "\n" +  "Ie_m1_random: " + str(Ie_m1_random)+  "\n" +  "w_n_nd_m_to_n_m1: " + str(w_n_nd_m_to_n_m1)+ "\n" +   "w_n_nd_c_to_n_pfc: " + str(w_n_nd_c_to_n_pfc) + "\n\n"
                                                                    
    time_run=mice.datetime.datetime.today()
    np.random.RandomState(seed=seed_mice)
    mice.scriviFile(initial_value, path, "initial_value_of_section.txt")
 
   
    for count in range(0,number_of_section):
    
        output_file_name = []
        
        loop_value = "wc: " + str(wc) + "\n" + "wc: " + str(wc) + "\n" + "wpfc: " + str(wpfc) + "\n" + "seed_mice: "+ str(seed_mice) + "\n" + "w_inc_pfc_i: "+ str(w_inc_pfc_i) + "\n" + "w_inc_pfc_i_nsp: "+ str(w_inc_pfc_i_nsp) + "\n" + "weight_k_increment_i_m: "+ str(weight_k_increment_i_m) + "\n" + "weight_k_increment_i_c: "+ str(weight_k_increment_i_c) + "\n" + "weight_gain_i: " + str(+weight_gain_i) + "\n" + "w_io_bm_i: "+str(w_io_bm_i) + "\n" + "w_io_bc_i: "+ str(w_io_bc_i) + "\n" + "w_cs_stim: "+ str(w_cs_stim)  + "\n" + "Ie_pfc_random: " + str(Ie_pfc_random) + "\n" + "Ie_m1_random: " + str(Ie_m1_random) + "\n" +  "w_n_nd_m_to_n_m1: " + str(w_n_nd_m_to_n_m1) + "\n" +  "w_n_nd_c_to_n_pfc: " + str(w_n_nd_c_to_n_pfc) +"\n\n"
        mice.scriviFile(loop_value, path, "value_of_section_loop.txt")
        
        out_net = mice.runNetwork (wm, wc, wpfc, seed_mice, w_inc_pfc_i, weight_gain_i, w_io_bm_i, w_io_bc_i, w_cs_stim, w_inc_pfc_i_nsp, weight_k_increment_i_c, weight_k_increment_i_m, Ie_pfc_random,Ie_m1_random, w_n_nd_m_to_n_m1, w_n_nd_c_to_n_pfc)


        # make the file name
        # file nr 0 granular cells line motor
        output_file_name.append( "output_gm_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 1 granular cells line cognitive
        output_file_name.append( "output_gc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 2 io line motor
        output_file_name.append( "output_im_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 3 io line cognitive
        output_file_name.append( "output_ic_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 4
        output_file_name.append( "output_m1_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 5
        output_file_name.append( "output_pfc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 6
        output_file_name.append( "output_pm_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 7
        output_file_name.append( "output_pc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 8 granular potential motor
        output_file_name.append( "output_vm_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 9 granular potential cognitive
        output_file_name.append( "output_vc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 10 granular spike motor
        output_file_name.append( "output_sm_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 11 granular spike cognitive
        output_file_name.append( "output_sc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 12
        output_file_name.append( "output_wm_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 13
        output_file_name.append( "output_wc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")        
        # file nr 14
        output_file_name.append( "output_wp_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 15
        output_file_name.append( "output_mc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #sender m1
        # file nr 16
        output_file_name.append( "output_px_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #sender pfc
        # file nr 17
        output_file_name.append( "output_pj_m"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #potential and time pj line motor       
        # file nr 18
        output_file_name.append( "output_pj_c"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #potential and time pj line cognitive  
        # file nr 19
        output_file_name.append( "output_ndm_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 20
        output_file_name.append( "output_ndc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")    

        
        # write  file
        for name_file in output_file_name: mice.scriviFile(name_file, path,"list_file_name.csv")
        #
        for gc in out_net["timesSpikeGC_m"]: mice.scriviFile(gc, path, output_file_name[0]) # output times of spike granular cells
        #
        for gc in out_net["timesSpikeGC_c"]: mice.scriviFile(gc, path, output_file_name[1]) # output times of spike granular cells
        #        
        for m1 in out_net["timesSpikeM1"]: mice.scriviFile(m1, path, output_file_name[4]) # output times of spike m1
        #    
        for pfc in out_net["timesSpikePFC"]: mice.scriviFile(pfc, path, output_file_name[5]) # output times of spike pfc
        #        
        for nd_m in out_net["timesSpikeND_m"]: mice.scriviFile(nd_m, path, output_file_name[19]) # output times of spike nd_m
        #    
        for nd_c in out_net["timesSpikeND_c"]: mice.scriviFile(nd_c, path, output_file_name[20]) # output times of spike nd_c
        #
        out_net_io_m=np.matrix([np.array(out_net["timesSpikeIO_m"]),np.array(out_net["sendersIO_m"])])
        np.savetxt(path+output_file_name[2], out_net_io_m.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_io_c=np.matrix([np.array(out_net["timesSpikeIO_c"]),np.array(out_net["sendersIO_c"])])
        np.savetxt(path+output_file_name[3], out_net_io_c.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_sgc_m=np.matrix([np.array(out_net["timesSpikeGC_m"]),np.array(out_net["sendersGC_m"])])
        np.savetxt(path+output_file_name[10], out_net_sgc_m.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_sgc_c=np.matrix([np.array(out_net["timesSpikeGC_c"]),np.array(out_net["sendersGC_c"])])
        np.savetxt(path+output_file_name[11], out_net_sgc_c.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_vgc_m=np.matrix([np.array(out_net["tPotentialGC_m"]),np.array(out_net["potentialGC_m"])])
        np.savetxt(path+output_file_name[8], out_net_vgc_m.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_vgc_c=np.matrix([np.array(out_net["tPotentialGC_c"]),np.array(out_net["potentialGC_c"])])
        np.savetxt(path+output_file_name[9], out_net_vgc_c.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_pj_m=np.matrix([np.array(out_net["timesSpikePj_m"]),np.array(out_net["sendersPj_m"])]) 
        np.savetxt(path+output_file_name[6], out_net_pj_m.transpose(), delimiter=",", fmt='%10.1f') 
        #
        out_net_pj_c=np.array([np.array(out_net["timesSpikePj_c"]),np.array(out_net["sendersPj_c"])])
        np.savetxt(path+output_file_name[7], out_net_pj_c.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_wm=np.array([np.array(out_net["tw_chart_motor"]),np.array(out_net["w_chart_motor"])])
        np.savetxt(path+output_file_name[12], out_net_wm.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_wc=np.array([np.array(out_net["tw_chart_cognitive"]),np.array(out_net["w_chart_cognitive"])])
        np.savetxt(path+output_file_name[13], out_net_wc.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_wpfc=np.array([np.array(out_net["t_chart_pfc"]),np.array(out_net["w_chart_pfc"])])
        np.savetxt(path+output_file_name[14], out_net_wpfc.transpose(), delimiter=",", fmt='%10.1f')

        out_net_m1c=np.array([np.array(out_net["timesSpikeM1"]),np.array(out_net["senderM1"])])
        np.savetxt(path+output_file_name[15], out_net_m1c.transpose(), delimiter=",", fmt='%10.1f')     
        
        out_net_pfc=np.array([np.array(out_net["timesSpikePFC"]),np.array(out_net["senderPFC"])])
        np.savetxt(path+output_file_name[16], out_net_pfc.transpose(), delimiter=",", fmt='%10.1f')
        
        out_net_pj2_m=np.matrix([np.array(out_net["tPotentialPj_m"]),np.array(out_net["potentialPj_m"])])
        np.savetxt(path+output_file_name[17], out_net_pj2_m.transpose(), delimiter=",", fmt='%10.1f')
        
        out_net_pj2_c=np.matrix([np.array(out_net["tPotentialPj_c"]),np.array(out_net["potentialPj_c"])])
        np.savetxt(path+output_file_name[18], out_net_pj2_c.transpose(), delimiter=",", fmt='%10.1f')
        
        
        
        wc=np.array(out_net["w_chart_cognitive"])[-1]
        wm=np.array(out_net["w_chart_motor"])[-1]
        wpfc=np.array(out_net["w_chart_pfc"])[-1]
        
