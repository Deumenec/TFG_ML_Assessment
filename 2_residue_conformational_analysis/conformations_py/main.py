########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Executing the conformation classification algorithm                                       #
#                                                                                                      #
########################################################################################################

import numpy as np
import os
import matplotlib.pyplot as plt
import MDAnalysis as mda 
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align

import useful_functions
import residue
import config
import __init__

aminoacid_ensemble = []

for aminoacid in config.aminoacids:
    aminoacid_ensemble.append(residue.Residue(aminoacid, config.len_aminoacids[aminoacid], config.source, config.conformations, config.output, config.threshold))


Train   = False     #Read into the conformational database the all the proteins stored in the dades folder
Read    = True      #Open the conformations saved in the conformacions folder
Save    = False     #Save the computed conformations into the conformacions folder
Analize = True     #Ask for molecule files to analize


if (Train == True):
    counter = 1
    for file in os.listdir(config.source+"/"):
        print(counter)
        counter += 1
        if file.endswith(".pdb"):
            u = mda.Universe(config.source+"/"+file)
            for i in range(len(aminoacid_ensemble)):
                aminoacid_ensemble[i].add_universe(u, os.path.splitext(file)[0],verbose = False) 
                
if (Read == True):
    for i in range(len(aminoacid_ensemble)):
        aminoacid_ensemble[i].open()
    
if (Save == True):
    for i in range(len(aminoacid_ensemble)):
        print(aminoacid_ensemble[i])
        aminoacid_ensemble[i].save()

if (Analize == True):
    while True:
        analized_prote  = input("quin fitxer de proteïna vols analitzar? Escriu break per parar: ")
        if analized_prote == "break":
            break
        u = mda.Universe(config.source_2+"/"+analized_prote)
        for i in range(len(aminoacid_ensemble)):
            freq_1, freq_2 = aminoacid_ensemble[i].frequencies_vector(u)
            useful_functions.ploter(freq_1, freq_2, config.output + "/" +os.path.splitext(analized_prote)[0], config.aminoacids[i]+os.path.splitext(analized_prote)[0], "pdf")

for i in range(len(aminoacid_ensemble)):
    aminoacid_ensemble[i].save()















