'''
@Author Piernadrea Mirino,  Daniele Caligiore,  Emiliano Trimarco  2020-2021
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
    
    wf = round(random.uniform(5., 5.), 2) #initializzation w kernel fast
    ws = round(random.uniform(20., 20.), 2) #initializzation w kernel slow 
    weight_k_increment_i_f=round(random.uniform(0.05, 0.05), 2) #0.05 costant w increment weight Granular cells to PJ no spike
    weight_k_increment_i_s=weight_k_increment_i_f

    w_cs_stim = round(random.uniform(500, 500), 2)#w CS stim 
    
    w_io_bf = random.randint(500, 500)*-1 #w US
    w_io_bf_i = w_io_bf #w io to pj  US
    w_io_bs_i = w_io_bf #w io to pj  US 
    
    weight_gain_i = round(random.uniform(5., 5.),2) #gain
    Ie_pfc_random = round(random.uniform(300.,365.),2)
    Ie_m1_random = round(random.uniform(300., 365.),2)
    
    w_n_nd_f_to_n_m1 = round(random.uniform(100., 100.), 2)
    w_n_nd_s_to_n_pfc = round(random.uniform(50., 50.), 2)
    
 
    path=os.getcwd()+"/session_mice_"+str(n_mice)+"/"
    
    mice.scriviFile(seed_mice, path, "seed.txt")
        
    ISI= 20
  
    
    #==============================================================================
    # 
    #==============================================================================
    
    
    initial_value = "wf: " + str(wf) + "\n" + "ws: " + str(ws) + "\n" + "wpfc: " + str(wpfc) + "\n" + "seed_mice: "+ str(seed_mice) + "\n" + "w_inc_pfc_i: "+ str(w_inc_pfc_i) + "\n" + "w_inc_pfc_i_nsp: "+ str(w_inc_pfc_i_nsp) + "\n" + "weight_k_increment_i_f: "+ str(weight_k_increment_i_f) + "\n" + "weight_k_increment_i_s: "+ str(weight_k_increment_i_s) + "\n" + "weight_gain_i: " + str(+weight_gain_i) + "\n" + "w_io_bf_i: "+str(w_io_bf_i) + "\n" + "w_io_bs_i: "+ str(w_io_bs_i) + "\n" + "w_cs_stim: "+ str(w_cs_stim) + "\n" + "Ie_pfc_random: " + str(Ie_pfc_random) + "\n" +  "Ie_m1_random: " + str(Ie_m1_random)+  "\n" +  "w_n_nd_f_to_n_m1: " + str(w_n_nd_f_to_n_m1)+ "\n" +   "w_n_nd_s_to_n_pfc: " + str(w_n_nd_s_to_n_pfc) + "\n\n"
                                                                    
    time_run=mice.datetime.datetime.today()
    np.random.RandomState(seed=seed_mice)
    mice.scriviFile(initial_value, path, "initial_value_of_section.txt")
 
   
    for count in range(0,number_of_section):
    
        output_file_name = []
        
        loop_value = "wfc: " + str(wf) + "\n" + "ws: " + str(ws) + "\n" + "wpfc: " + str(wpfc) + "\n" + "seed_mice: "+ str(seed_mice) + "\n" + "w_inc_pfc_i: "+ str(w_inc_pfc_i) + "\n" + "w_inc_pfc_i_nsp: "+ str(w_inc_pfc_i_nsp) + "\n" + "weight_k_increment_i_f: "+ str(weight_k_increment_i_f) + "\n" + "weight_k_increment_i_s: "+ str(weight_k_increment_i_s) + "\n" + "weight_gain_i: " + str(+weight_gain_i) + "\n" + "w_io_bf_i: "+str(w_io_bf_i) + "\n" + "w_io_bs_i: "+ str(w_io_bs_i) + "\n" + "w_cs_stim: "+ str(w_cs_stim)  + "\n" + "Ie_pfc_random: " + str(Ie_pfc_random) + "\n" + "Ie_m1_random: " + str(Ie_m1_random) + "\n" +  "w_n_nd_f_to_n_m1: " + str(w_n_nd_f_to_n_m1) + "\n" +  "w_n_nd_s_to_n_pfc: " + str(w_n_nd_s_to_n_pfc) +"\n\n"
        mice.scriviFile(loop_value, path, "value_of_section_loop.txt")
        
        out_net = mice.runNetwork (wf, ws, wpfc, seed_mice, w_inc_pfc_i, weight_gain_i, w_io_bf_i, w_io_bs_i, w_cs_stim, w_inc_pfc_i_nsp, weight_k_increment_i_s, weight_k_increment_i_f, Ie_pfc_random,Ie_m1_random, w_n_nd_f_to_n_m1, w_n_nd_s_to_n_pfc)


        # make the file name
        # file nr 0 granular cells line fast
        output_file_name.append( "output_gf_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 1 granular cells line slow
        output_file_name.append( "output_gs_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 2 io line fast
        output_file_name.append( "output_if_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 3 io line slow
        output_file_name.append( "output_is_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 4
        output_file_name.append( "output_m1_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 5
        output_file_name.append( "output_pfc_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 6
        output_file_name.append( "output_p1_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 7
        output_file_name.append( "output_p2_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 8 granular potential fast
        output_file_name.append( "output_vf_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 9 granular potential slow
        output_file_name.append( "output_vs_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 10 granular spike fast
        output_file_name.append( "output_sf_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 11 granular spike slow
        output_file_name.append( "output_ss_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 12
        output_file_name.append( "output_wf_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 13
        output_file_name.append( "output_ws_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")        
        # file nr 14
        output_file_name.append( "output_wp_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")
        # file nr 15
        output_file_name.append( "output_ms_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #sender m1
        # file nr 16
        output_file_name.append( "output_ps_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #sender pfc
        # file nr 17
        output_file_name.append( "output_pA_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #potential and time pj line fast       
        # file nr 18
        output_file_name.append( "output_pL_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv") #potential and time pj line slow  
        # file nr 19
        output_file_name.append( "output_ndf_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")############
        # file nr 20
        output_file_name.append( "output_nds_"+ str(time_run) +"_ISI_"+str(ISI)+"Session_"+str(count)+".csv")###############     

        
        # write  file
        for name_file in output_file_name: mice.scriviFile(name_file, path,"list_file_name.csv")
        #
        for gc in out_net["timesSpikeGC_f"]: mice.scriviFile(gc, path, output_file_name[0]) # output times of spike granular cells
        #
        for gc in out_net["timesSpikeGC_s"]: mice.scriviFile(gc, path, output_file_name[1]) # output times of spike granular cells
        #        
        for m1 in out_net["timesSpikeM1"]: mice.scriviFile(m1, path, output_file_name[4]) # output times of spike m1
        #    
        for pfc in out_net["timesSpikePFC"]: mice.scriviFile(pfc, path, output_file_name[5]) # output times of spike pfc
        #        
        for ndf in out_net["timesSpikeNDf"]: mice.scriviFile(ndf, path, output_file_name[19]) # output times of spike ndf##################
        #    
        for nds in out_net["timesSpikeNDs"]: mice.scriviFile(nds, path, output_file_name[20]) # output times of spike nds################
        #
        out_net_io_f=np.matrix([np.array(out_net["timesSpikeIO_f"]),np.array(out_net["sendersIO_f"])])
        np.savetxt(path+output_file_name[2], out_net_io_f.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_io_s=np.matrix([np.array(out_net["timesSpikeIO_s"]),np.array(out_net["sendersIO_s"])])
        np.savetxt(path+output_file_name[3], out_net_io_s.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_sgc_f=np.matrix([np.array(out_net["timesSpikeGC_f"]),np.array(out_net["sendersGC_f"])])
        np.savetxt(path+output_file_name[10], out_net_sgc_f.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_sgc_s=np.matrix([np.array(out_net["timesSpikeGC_s"]),np.array(out_net["sendersGC_s"])])
        np.savetxt(path+output_file_name[11], out_net_sgc_f.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_vgc_f=np.matrix([np.array(out_net["tPotentialGC_f"]),np.array(out_net["potentialGC_f"])])
        np.savetxt(path+output_file_name[8], out_net_vgc_f.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_vgc_s=np.matrix([np.array(out_net["tPotentialGC_s"]),np.array(out_net["potentialGC_s"])])
        np.savetxt(path+output_file_name[9], out_net_vgc_s.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_pjf=np.matrix([np.array(out_net["timesSpikePjF"]),np.array(out_net["sendersPjF"])]) #p1
        np.savetxt(path+output_file_name[6], out_net_pjf.transpose(), delimiter=",", fmt='%10.1f') 
        #
        out_net_pjs=np.array([np.array(out_net["timesSpikePjS"]),np.array(out_net["sendersPjS"])])#p2
        np.savetxt(path+output_file_name[7], out_net_pjs.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_wf=np.array([np.array(out_net["tw_chart_fast"]),np.array(out_net["w_chart_fast"])])
        np.savetxt(path+output_file_name[12], out_net_wf.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_ws=np.array([np.array(out_net["tw_chart_slow"]),np.array(out_net["w_chart_slow"])])
        np.savetxt(path+output_file_name[13], out_net_ws.transpose(), delimiter=",", fmt='%10.1f')
        #
        out_net_wpfc=np.array([np.array(out_net["t_chart_pfc"]),np.array(out_net["w_chart_pfc"])])
        np.savetxt(path+output_file_name[14], out_net_wpfc.transpose(), delimiter=",", fmt='%10.1f')

        out_net_m1s=np.array([np.array(out_net["timesSpikeM1"]),np.array(out_net["senderM1"])])
        np.savetxt(path+output_file_name[15], out_net_m1s.transpose(), delimiter=",", fmt='%10.1f')     
        
        out_net_pfc=np.array([np.array(out_net["timesSpikePFC"]),np.array(out_net["senderPFC"])])
        np.savetxt(path+output_file_name[16], out_net_pfc.transpose(), delimiter=",", fmt='%10.1f')
        
        out_net_PjA=np.matrix([np.array(out_net["tPotentialPjA"]),np.array(out_net["potentialPjA"])])
        np.savetxt(path+output_file_name[17], out_net_PjA.transpose(), delimiter=",", fmt='%10.1f')
        
        out_net_PjL=np.matrix([np.array(out_net["tPotentialPjL"]),np.array(out_net["potentialPjL"])])
        np.savetxt(path+output_file_name[18], out_net_PjL.transpose(), delimiter=",", fmt='%10.1f')
        
        
        
        ws=np.array(out_net["w_chart_slow"])[-1]
        wf=np.array(out_net["w_chart_fast"])[-1]
        wpfc=np.array(out_net["w_chart_pfc"])[-1]
        
