########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Defining the Residue class                                                                #
#                                                                                                      #
########################################################################################################

import config
import MDAnalysis as mda 
from MDAnalysis.analysis import align
import numpy as np
import os

class Residue:
    def __init__(self, residue_id, residue_len, source_dir, conformations_dir, output_dir, threshold):
        self.residue_id          = residue_id
        self.residue_len         = residue_len
        self.source_dir          = source_dir
        self.conformations_dir   = conformations_dir
        self.output_dir          = output_dir
        self.threshold           = threshold*residue_len**0.5
        self.conformations       = []
        self.total_num           = 0             #Number of diferent conformations found
        self.residue_num         = 0             #Total of conformations analyzed
        
    class Conformation:
        def __init__(self, conformation_id, conformation_cord, origin, instances_num = 1):
            self.conformation_id  = conformation_id
            self.conformation_cord = conformation_cord
            self.instances_num = instances_num
            self.instances = [conformation_cord] #Aquí podria guardar totes les conformacions observades amb l'RMSD fixat però va massa lent
            self.origin = origin #Saves the first file where the conformation was found...could save all the files where it was found!!!
                        
        def compare(self, compare_conformation, compare_threshold):
            "Check if a conformation is equal within the given RMSD value"
            if(len(compare_conformation.atoms)!= len(self.conformation_cord.atoms)):
                #print("Hi ha una "+self.conformation_id+" de len "+str(len(compare_conformation))+" en comptes de " + str(len(self.conformation_cord))) Per excloure les conformacions terminals i d'altres errors
                return
            _, aligned_rmsd = align.alignto(compare_conformation, self.conformation_cord) 
            if aligned_rmsd <= compare_threshold:
                self.instances_num +=1
                return True
            return False
    
        def __repr__(self):
            return f"Conformation of '{self.conformation_id}' with '{self.instances_num}' instances first seen in '{self.origin}'" 
        
    def check_conformation(self, checking_conformation, check_name):
        "Given a conformation, it checks if it coincides with an already stored conformation. If that is the case, it stores it within that bin. Otherwise, it creates a conformation class for that conformation."
        if (len(checking_conformation.atoms) != self.residue_len):
            return
        for check_stored_conformation in self.conformations:
            if check_stored_conformation.compare(checking_conformation, self.threshold) == True:
                check_stored_conformation.instances.append(checking_conformation)
                return
        self.conformations.append(self.Conformation(self.residue_id, checking_conformation, check_name))
        #print("halou"+self.residue_id)
        return

    def add_universe(self, universe, name, verbose = False):
        "Checks the conformation for all instances of a residue inside a universe"
        add_all_residues= universe.select_atoms("protein and not (name OXT or name NTER) and resname "+self.residue_id)
        add_frames = add_all_residues.residues
        for add_frame in add_frames:
            if len(self.conformations)==0:
                if len(add_frame.atoms)== self.residue_len:
                    self.total_num +=1
                    self.conformations.append(self.Conformation(self.residue_id, add_frame, name))
            else:
                self.total_num +=1
                self.check_conformation(add_frame, name)
        if verbose == True:
            print(self)
        return
    def frequencies_vector(self, freq_universe):
        "Given a universe, it adds it and returns aminoacid frequencies, printing the performed calculation "
        print(self.residue_id)
        freq_frequencies_0 = []
        freq_frequencies_1 = []
        for freq_conformation in self.conformations:
            freq_frequencies_0.append(freq_conformation.instances_num)
        self.add_universe(freq_universe, "analysis") #Afegeix el nou univers que es vol analitzar
        for freq_conformation in self.conformations:
            freq_frequencies_1.append(freq_conformation.instances_num)
        return freq_frequencies_0, freq_frequencies_1
            
    def open(self):
        "Opens conformations saved in conformations_dir, no és necessari realment així que no ho he programat aquí"
        open_conformations = [f for f in os.listdir(self.conformations_dir +"/"+self.residue_id ) if not f.endswith(".txt")]
        
        #Resets the open conformations
        self.conformations       = []
        self.total_num           = 0             
        self.residue_num         = 0
        
        for open_conformation in range(len(open_conformations)):
            with open (self.conformations_dir +"/"+self.residue_id +"/"+str(open_conformation)+"data.txt") as file:
                lines = file.read()
                lines = lines.splitlines()
            open_coordinates = mda.Universe(self.conformations_dir +"/"+self.residue_id +"/"+str(open_conformation)+".pdb")
            self.conformations.append(self.Conformation(self.residue_id, open_coordinates, lines[1], instances_num = int(lines[0])))
            del open_coordinates #Per esborrar l'univers un cop s'han guardat les dades
    
    def save(self):
        "Saves the conformations in conformations_dir"
        os.makedirs(self.conformations_dir + "/" + self.residue_id, exist_ok=True)
        for save_conformation in range(len(self.conformations)):
            with mda.Writer(self.conformations_dir + "/" + self.residue_id + "/"+str(save_conformation)+'.pdb', self.residue_len) as file:
                save_conf= self.conformations[save_conformation].conformation_cord
                save_conf.atoms.translate(-save_conf.atoms.center_of_mass())
                file.write(save_conf)
            with open(self.conformations_dir + "/" + self.residue_id + "/"+str(save_conformation)+'data.txt', "w") as file:
                file.write(str(self.conformations[save_conformation].instances_num)+ "\n" +self.conformations[save_conformation].origin)
        return

    def __repr__(self):
        return f"Aminoacid: \nName='{self.residue_id}' \nNumber of atoms in residue= '{self.residue_len}' \nNumber of conformations found='{len(self.conformations)}'\nTotal of conformations analyzed='{self.total_num}'\nRMSD Threshhold used= '{self.threshold}'\n "
        