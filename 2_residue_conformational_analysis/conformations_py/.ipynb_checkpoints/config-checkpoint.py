########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Defining the main parameters for the method used and file operations                      #
#                                                                                                      #
########################################################################################################

#Allowed RMSD treshold for residues in a given bin.
threshold = 0.2  

#Directory for reading data
source = "dades"

#Directory for storing conformations
conformations = "conformacions"

#Directory for outputing data
output = "resultats"

#Names of the aminoacids present
aminoacids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", 
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]

#number of atoms expected in each residue (to aviod conflicts with terminal residues)
len_aminoacids = {
"ALA" : 5, "ARG" : 11, "ASN" : 8, "ASP" : 8, "CYS" : 6, "GLN" : 9, "GLU" : 9, "GLY" : 4, "HIS" : 10, "ILE" : 8, "LEU" : 8, "LYS" : 9, "MET" : 8, "PHE" : 11, "PRO" : 7, "SER" : 6, "THR" : 7, "TRP" : 14, "TYR" : 12, "VAL" : 7
}