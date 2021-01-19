#*************************************************************
#*************************************************************
#*************************************************************
#*************************************************************
#*************************************************************
#*************************************************************
import numpy as np
import pylab as pl
import sys
import csv
from collections import Counter
import random


n_subject = int(sys.argv[3])
n_session = 11
peak_array=np.ones((n_subject,n_session))


ft = 1
inz= 0
fin= 550

        
pl.rcParams.update({'figure.max_open_warning': 0})

try:
    base_path=sys.argv[1]
    session_number=int(sys.argv[3])
    #print sys.argv
except:
    print ("Manca il path della cartella che contiene le sessioni o il numero delle sessioni \n")
    print ("Sample: python <name of fil,e:.py> '</directory path/:str>' <number of session start:int> <number of session stop:int>")
    exit()


def peak(b, ft):
	quotients = [number / ft for number in b]
	#print(quotients)
	cc = [ round(num) for num in quotients ]
	count = Counter(cc)
	freq_list = count.values()

	max_cnt = max(freq_list)
	total = freq_list.count(max_cnt)

	most_common = count.most_common(total)


	a = [elem[0] for elem in most_common]
	print(np.dot(a,ft))	
	s = sum(a) / len(a)
	res = round(s*ft)
	return res     


def make_chart (path, mice=0, inz=0, fin=550): 
          
    path=base_path+path
    print ("\n\n")    
    print ("***/%/*** Run session => " + path + "  ***/%/***")
    print ("\n\n")
                                       
    with open(path+'list_file_name.csv') as csv_name_file:
        csv_reader_name_file = csv.reader(csv_name_file, delimiter=',')
        data=[]
        counter = []
        CC = []
        el = []
        array_massimo = []
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
                
                print ("Import file <= " + path+row_name_file[0] + "\n")
                
                el = [x for x in data if  inz < x < fin] 

                #print(el)
                
                try:
                    CC = peak(el,ft)
                    array_massimo.append(CC)
                    print(CC)
                    CC= []
                except:
                    array_massimo.append(round(random.uniform(0, 0), 2))
                    print("0000000000000000000000000000000000____ERROR____00000000000000000000000000000000000")
                    
                                                         
                number_session.append(int(str.replace(row_name_file[0][-6:-4],'_','')))
          
            row=[]
            data=[]
            line_count_name_file += 1
    
    print(array_massimo)
    peak_array[i]= array_massimo
    array_massimo=[]
    #print(peak_array)
    return peak_array






    
for i in range (int(sys.argv[2]), int(sys.argv[3])):
	
	
	peak_arrayCC = make_chart('control/session_mice_'+str(i)+'/', mice=i,inz=0, fin=550)                

peak_array=np.ones((n_subject,n_session))


for i in range (int(sys.argv[2]), int(sys.argv[3])):
	
	
	peak_arraySS = make_chart('ASD/session_mice_'+str(i)+'/', mice=i, inz=0, fin= 500)

peak_array=np.ones((n_subject,n_session))	


peak_arrayC = np.delete(peak_arrayCC, 0, 1)

peak_mean = np.mean(peak_arrayC, axis= 0)
stddev_peak_array = np.std(peak_arrayC, axis= 0)

print("==========PEAK value===============")
print(peak_arrayC)
#print("=================================")
#print("==========PEAK mean===============")
#print(peak_mean)
#print("=================================")
#print("==========STD value===============")
#print(stddev_peak_array)
#print("=================================")


peak_arrayS = np.delete(peak_arraySS, 0, 1)
peak_meanS = np.mean(peak_arrayS, axis= 0)
stddev_peak_arrayS = np.std(peak_arrayS, axis= 0)

print("==========PEAK valueSSS===============")
print(peak_arrayS)
"""
print("=================================")
print("==========PEAK meanSSS===============")
print(peak_meanS)
print("=================================")
print("==========STD valueSSS===============")
print(stddev_peak_arrayS)
print("=================================")


"""

        
x = np.arange(1,(n_session), 1)
    
pl.errorbar(x, peak_mean, stddev_peak_array, marker='o', markerfacecolor='white', markersize=4, color='grey', linewidth=1, label='Control')
pl.errorbar(x, peak_meanS, stddev_peak_arrayS, marker='o', markerfacecolor='black', markersize=4, color='grey', linewidth=1, label='ASD')

pl.suptitle('Peak Eyeblink Latency', fontsize=16)
pl.xlabel('Session')
pl.ylabel('Peak Latency (msec)')
pl.legend(loc = 1)
pl.yticks(np.arange(250, 475, 25))
pl.xticks(x)
pl.show()
    


