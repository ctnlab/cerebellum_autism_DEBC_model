# libraries

import numpy as np
import pylab as pl
import sys
import matplotlib.pyplot as plt


n_subject = 15
n_session = 11
control= np.ones((n_subject,n_session))
sperimental= np.ones((n_subject,n_session))


def carica (path):
    base_path=sys.argv[1]
    countRow = 0
    txt_reader_row =[]
    data = []
    path=base_path+path
    file_caricati = np.ones((n_subject,n_session))


    for h in range(0, n_subject):
        
        txt_reader = []
        with open(path+str(h)+".txt") as txt_file:
           
            txt_reader_row = np.loadtxt(txt_file, delimiter=',')
            data = []
            for row in txt_reader_row:

                #print(countRow)
                
                data.append(float(row))
                countRow=countRow+1
                print (data) 
            file_caricati[h] = data                           
            
    return file_caricati
    
for i in range (int(sys.argv[2]), int(sys.argv[3])): 
   
    
    controlC = carica('control/CR-CS_US/cr_medio_di_sessione_ession_mice_') 
    
controlT= np.mean(control, axis=1)    
print("=========================control=======TTT======")
print(controlT)    
    
for i in range (int(sys.argv[2]), int(sys.argv[3])): 
   
    
    sperimentalS = carica('ASD/CR-CS_US/cr_medio_di_sessione_ession_mice_') 

control= np.delete(controlC, 0, 1)
sperimental=np.delete(sperimentalS, 0, 1)


control1 = np.mean(control, axis=0)
sd_control = np.std(control, axis= 0) 

sperimental1 = np.mean(sperimental, axis=0)
sd_sperimental=np.std(sperimental, axis= 0) 

controlT= np.mean(control, axis=1)
sperimentalT= np.mean(sperimental, axis=1)

print(control)
print (sperimental)
print(sd_control)
print (sd_sperimental)
x = np.arange(1,(n_session), 1)


 
# multiple line plot
pl.errorbar( x, control1, sd_control, marker='o', markerfacecolor='white', markersize=4, color='grey', linewidth=1, label='Control')
pl.errorbar( x, sperimental1, sd_sperimental, marker='o', markerfacecolor='black', markersize=4, color='grey', linewidth=1, label='Autism')
pl.suptitle('Percent CRs', fontsize=16)
pl.xlabel('Session')
pl.ylabel('Percent CR')
pl.yticks(np.arange(0, 120, 20))
pl.xticks(x)
pl.legend(loc = 4)




print("=========================control=======TTT======")
print(control1)
print ("=================sperimental==========TTT======")
print (sperimental1)

pl.show() 
