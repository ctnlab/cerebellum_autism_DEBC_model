import pylab as pl
import elephant as el
import numpy as np
import math as m
import time
import sys
import datetime
import quantities as pq
import neo
import scipy
from matplotlib.backends.backend_pdf import PdfPages
from time import sleep
import csv
from os.path import exists
import os

#firing rate

def firingRate(spike_train,sig=100,  start=0., stop=400., ch=2.):
    
    kernel_gau=el.kernels.GaussianKernel(sigma=sig*pq.ms)
    
    if ch == 1.: #all
        firing_rate=el.statistics.instantaneous_rate(neo.SpikeTrain(spike_train, units=pq.ms, t_start=np.array(spike_train).min(), t_stop=np.array(spike_train).max()), sampling_period=pq.Quantity([1], 'ms'), kernel=kernel_gau)
  
    if ch == 2.: #period
          
        SpikeTrainL = list(np.where(np.logical_and(np.array(spike_train)>=start, np.array(spike_train)<=stop)))     
        array_a =neo.SpikeTrain ( SpikeTrainL, units=pq.ms, t_start=np.array(SpikeTrainL).min(),  t_stop=np.array(SpikeTrainL).max())
        firing_rate=el.statistics.instantaneous_rate(array_a, sampling_period=pq.Quantity([1], 'ms'), kernel=kernel_gau)
        
     
    return firing_rate # Write  file
     

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


def calculate_cr (path, mice = 0):
    
    version = 8.2
    label_grf=''
    session_time = 2500
    FRm_array=[]
    CRp_array=[]
    FR_mean_range_of_value_m1=[]
    FR_mean_range_of_value_m1_after_first_cs = []
    
    try:
        base_path=sys.argv[1]
        session_number=sys.argv[2]
        #print sys.argv
    except:
        print "The path of the folder containing the sessions or the number of sessions is missing \n"
        print "Sample: python <name of file:.py> '</directory path/:str>' <number of session:int>"
        exit()
  
    path=base_path+path
    
    print "\n\n"    
    print "***/%/*** Run session => " + path + "  ***/%/***"
    print "\n\n"
                            
    #make a python dictionary with firingRate of value find into file name content'list_file_name.csv'                    
    with open(path+'list_file_name.csv') as csv_name_file:
        csv_reader_name_file = csv.reader(csv_name_file, delimiter=',')
        data=[]
        n_spike=[]
        fr_max = []
        firing_rate_m1_learning=[]
        learning=[]
        firing_rate_area= {}
        spike_area={}
        line_count_name_file = 0
        number_session=[]
        
        for row_name_file in csv_reader_name_file:
            

            try:
       
                with open(path+row_name_file[0]) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    line_count = 0
                 
                    for row in csv_reader:
                        data.append(float(row[0]))
                        line_count += 0                
            except:
                continue
            
     
            if (row_name_file[0][:9] == 'output_m1'): 
                
                out_fr=np.array(firingRate(data, ch=2.)).max()
                
                print (out_fr)
                
                firing_rate_m1_learning.append(out_fr)



            if (row_name_file[0][:9]=='output_m1' or row_name_file[0][:9]=='output_io' or row_name_file[0][:9]=='output_sg'):
                                               
                print ("Import file <= " + path+row_name_file[0] + "\n")                   
                firing_rate_area[row_name_file[0][:9]+"_"+str.replace(row_name_file[0][-6:-4],'_','')] = firingRate(data) # put firing rate into dictionary                       
                number_session.append(int(str.replace(row_name_file[0][-6:-4],'_','')))                    
                try:
                    spike_area[row_name_file[0][:9]+"_"+str.replace(row_name_file[0][-6:-4], '_', '')] = np.loadtxt(path+row_name_file[0], delimiter=',')
                except:
                    continue
                   
            row=[]
            data=[]
            line_count_name_file += 1


    with open(base_path+"fr_max.txt") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
     
        for row in csv_reader:
            fr_max.append(float(row[0]))

                                        
    for count_us in np.unique(number_session):
    
            
            print ("\n\n ========================================================================= \n")
            print ("\n\n ========================================================================= \n\n")              
            print ("====> SECTION NUMBER ---> " + str(count_us) + "\n\n")            

            
    
            CR = firing_rate_m1_learning[count_us]/np.array(fr_max).max()*100                            

            scriviFile(round(CR), base_path+"CR-CS_US/", "cr_medio_di_sessione_"+str(path[-14:-1])+".txt")
            
            CRp_array.append(round(CR))
            

            
                        
for i in range (0,int(sys.argv[2])):
    
    print i
        
    calculate_cr ('session_mice_'+str(i)+'/', mice=i) 
