import MDAnalysis as mda
from MDAnalysis.analysis import align
import warnings
import numpy as np
import matplotlib.pyplot as plt
import os

warnings.filterwarnings('ignore') # suppress some MDAnalysis warnings about PSF files
print("Using MDAnalysis version", mda.__version__)


aminoacids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", 
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]
results_folder = "resultats"
data_folder = "dades"
molecule_file = input("quin fitxer de molecules vols analitzar? ")
molecule_name = molecule_file.split(".")[0]

def folder_check(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

folder_check(results_folder)
folder_check(results_folder + "/" + molecule_name)

#funcions definides

def average_structure(av_residu, av_univers):
    """
    For a given aminoacid and a universe, returns the average shape of that aminoacid in the universe
    """
    av_univers = av_univers.select_atoms("resname "+av_residu)
    n_frames = len(av_univers.residues)
    n_atoms = len(av_univers.residues[0].atoms)
    with mda.Writer(results_folder + "/" + molecule_name+ "/" +av_residu +'.xtc', n_atoms) as w:
        for ts in range(n_frames):
            w.write(av_univers.residues[ts].atoms)
    av_univers=mda.Universe(results_folder + "/" + molecule_name+ "/"+av_residu +'.xtc')
    average = align.AverageStructure(av_univers,
                                     ref_frame=0).run()
    ref = average.results.universe
    ref.atoms.write(results_folder + "/" + molecule_name+ "/" + av_residu +"average.pdb")

def distortion_distribution(dist_univers, dist_residu, dist_reference):
    """
    for a given aminoacid it calclates the RMSD 
    """
    dist_univers = dist_univers.select_atoms("resname "+dist_residu)
    distribution= []
    for res_i in dist_univers.residues:
        print(dist_reference, res_i.atoms)
        distribution.append(mda.analysis.rms.RMSD(dist_reference, res_i.atoms))
    print(distribution)
def max_list(m_list, m_select):
    """
    Given a list of tuples it returns the one that has the max value at a given position
    """
    m_max= m_list[0]
    for element in m_list:
        if element[m_select]>m_max[m_select]:
            m_max = element
    return m_max

def n_alpha_dist(protein):
    """
    For a given protein, it returns the distance between the c_alpha atom and its -COO atached carbon for each residue
    """
    if not isinstance(protein, mda.core.universe.Universe):
        raise TypeError(f"Expected input_value to be of type mda.core.universe.Universe, but got {type(protein).__name__}")
    distances = []
    for res in protein.residues:
        nitrogen = res.atoms.select_atoms("name N")
        carbon_a = res.atoms.select_atoms("name CA")
        if len(nitrogen) == 1 and len(carbon_a) == 1:
            distance = np.linalg.norm(nitrogen.positions[0] - carbon_a.positions[0])
            distances.append((res.resid, res.resname, distance))
        else:
            print(f"Skipping residue {res.resid} ({res.resname}) due to lacking infotmation")
    return distances


###Main

u = mda.Universe(data_folder+"/"+molecule_file)
#Es determinen les estructures mitjanes
for aminoacid in aminoacids:
    average_structure(aminoacid, u)
#Es determinen les RMSD respecte les estructures mitjanes
for aminoacid in aminoacids:
    distortion_distribution(u, aminoacid, mda.Universe(results_folder + "/" + molecule_name+ "/" + aminoacid +"average.pdb").atoms)


    