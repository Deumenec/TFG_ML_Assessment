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
    def __init__(self, residue_id, source_dir, conformations_dir, output_dir, threshold):
        self.residue_id          = residue_id
        self.source_dir          = source_dir
        self.conformations_dir   = conformations_dir
        self.output_dir          = output_dir
        self.threshold           = threshold
        self.conformations       = []
        self.conformations_num   = 0
        
    class Conformation:
        def __init__(self, conformation_id, conformation_cord):
            self.conformation_id  = conformation_id
            self.conformation_cord = conformation_cord
            self.instances = [conformation_id]
            
        def compare(self, compare_conformation, compare_threshhold):
            "Check if a conformation is equal within the given RMSD value"
            _, aligned_rmsd = align.alignto(compare_conformation, self.conformation_cord) 
            if aligned_rmsd <= compare_threshold:
                return True
            return
        

    def check_conformation(checking_conformation):
        "Given a conformation, it checks if it coincides with an already stored conformation. If that is the case, it stores it within that bin. Otherwise, it creates a conformation class for that conformation."
        for check_stored_conformation in self.conformations:
            if check_stored_conformation.compare(checking_conformation, threshold) == True:
                check_stored_conformation.instances.append(checking_conformation)
                return
        self.conformations.append(Conformation(self.residue_id, checking_conformation))
        self.conformations_num += 1    
        return

    
    def add_universe(self, residue_source, data_directory, universe):
        "Checks the conformation for all instances of a residue inside a universe"
        add_all_residues= universe.select_atoms("protein and not (name OXT or name NTER) and resname"+self.residue_id)
        add_frames = add_all_residues.residues
        if len(add_frames) ==0:
            return
        for add_frame in add_frames:
            check_conformation(add_frame.atoms)

            
        
    
    def __repr__(self):
        return f"Aminoacid: \nName='{self.residue_id}' \nNumber of conformations='{len(self.conformations)}')"
        