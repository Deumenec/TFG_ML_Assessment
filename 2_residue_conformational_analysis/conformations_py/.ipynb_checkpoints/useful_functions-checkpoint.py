########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Definig useful functions used generally in all the program like plotting graphs           #
#                                                                                                      #
########################################################################################################

import os

aminoacids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", 
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]

def folder_check(folder_name):
    """
    Checks for a required folder required in the code ans
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
