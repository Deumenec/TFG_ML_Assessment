# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 18:26:01 2024

@author: dhe02
"""

import MDAnalysis as mda

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

from matplotlib import pyplot as plt

print("Using MDAnalysis version", mda.__version__)

structure = mda.Universe("1fqy.pdb")
bonds = structure.atoms.guess_bonds()
print(bonds)
alphas = structure.select_atoms("name CA")
names = structure.residues.resnames
