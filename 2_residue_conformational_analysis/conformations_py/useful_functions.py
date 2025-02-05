########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Definig useful functions used generally in all the program like plotting graphs           #
#                                                                                                      #
########################################################################################################

import os

def folder_check(folder_name):
    """
    Checks for a required folder required in the code ans
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

def add_file(prote):
    "Adds a file to the database of conformations"
    a