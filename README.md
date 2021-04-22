# cerebellum_autism_DEBC_model
Cortico-Cerebellar Hyper-Connections and Reduced Purkinje Cells Behind Abnormal Eyeblink Conditioning in Autism

To run the simulation, run the experiment file of each group. The software will generate activity files for all areas.

After generating the experiment file, within the folders of the two groups create a file "fr_max.txt" that contains the ideal firing rate values of M1 with respect to 100% correct responses during all sessions. Create this file before launching the CR calculation file. Use the cr_calcolo_control.py file to calculate the CR of the control group and cr_calcolo_ASD.py, to calculate the CR of the ASD group.
Sample: python <name of file: .py> ' < / directory path/ :str> ' < number of session: int>

When you start the group CR calculation, a "CR-CS_US" folder will be created inside the relative group directory.
Use the file "cr2.py" to make the graph of the two groups' CR.
Sample: python < name of file: .py> '< / directory path / :str >' < number of session start :int> < number of session stop :int>

Use the file "peak_value2.py" to make the peak latency plot of the two groups' responses.
Sample: python < name of file: .py> '< / directory path / :str >' < number of session start :int> < number of session stop :int>
