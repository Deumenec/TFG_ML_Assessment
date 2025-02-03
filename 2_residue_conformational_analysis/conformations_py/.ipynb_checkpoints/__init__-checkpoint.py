########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Initializes the working directory creating required folders                               #
#                                                                                                      #
########################################################################################################

import os
import config

#Set the working directory to the correct folder
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(script_dir, ".."))
os.chdir(parent_dir)

os.makedirs(config.source, exist_ok=True)
os.makedirs(config.output, exist_ok=True) 
os.makedirs(config.conformations, exist_ok=True)

