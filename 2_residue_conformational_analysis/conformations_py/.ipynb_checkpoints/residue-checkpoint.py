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

class Residue:
    def __init__(self, residue_id, residue_len, source_dir, conformations_dir, output_dir, threshold):
        self.residue_id          = residue_id
        self.residue_len         = residue_len
        self.source_dir          = source_dir
        self.conformations_dir   = conformations_dir
        self.output_dir          = output_dir
        self.threshold           = threshold
        self.conformations       = []
        self.total_num   = 0             #Number of diferent conformations found
        self.residue_num         = 0             #Total of conformations analyzed
        
    class Conformation:
        def __init__(self, conformation_id, conformation_cord):
            self.conformation_id  = conformation_id
            self.conformation_cord = conformation_cord
            self.instances_num = 1
            self.instances = [conformation_cord]
            
            
        def compare(self, compare_conformation, compare_threshold):
            "Check if a conformation is equal within the given RMSD value"
            if(len(compare_conformation)!= len(self.conformation_cord)):
                #print("Hi ha una "+self.conformation_id+" de len "+str(len(compare_conformation))+" en comptes de " + str(len(self.conformation_cord)))
                return
            _, aligned_rmsd = align.alignto(compare_conformation, self.conformation_cord) 
            if aligned_rmsd <= compare_threshold:
                self.instances_num +=1
                return True
            return False
        def __repr__(self):
            return f"Conformation of {self.conformation_id} with {self.instances_num} instances" 
        

    def check_conformation(self, checking_conformation):
        "Given a conformation, it checks if it coincides with an already stored conformation. If that is the case, it stores it within that bin. Otherwise, it creates a conformation class for that conformation."
        for check_stored_conformation in self.conformations:
            if check_stored_conformation.compare(checking_conformation, self.threshold) == True:
                check_stored_conformation.instances.append(checking_conformation)
                return
        self.conformations.append(self.Conformation(self.residue_id, checking_conformation))
        #print("halou"+self.residue_id)
        return

    
    def add_universe(self, universe):
        "Checks the conformation for all instances of a residue inside a universe"
        add_all_residues= universe.select_atoms("protein and not (name OXT or name NTER) and resname "+self.residue_id)
        add_frames = add_all_residues.residues
        for add_frame in add_frames:
            if len(self.conformations)==0:
                if len(add_frame.atoms)== self.residue_len:
                    self.total_num +=1
                    self.conformations.append(self.Conformation(self.residue_id, add_frame.atoms))
            else:
                self.total_num +=1
                self.check_conformation(add_frame.atoms)
        return
     
    def __repr__(self):
        return f"Aminoacid: \nName='{self.residue_id}' \nNumber of atoms in residue= '{self.residue_len}' \nNumber of conformations found='{len(self.conformations)}'\nTotal of conformations analyzed='{self.total_num}'\n"
        