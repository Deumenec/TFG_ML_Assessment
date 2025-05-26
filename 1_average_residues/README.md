3/2/2025
This is the first program wrote in the TFG with the goal of validating aminoacid structure proportioned by AF3. It calculates the average shape for each aminoacid by aligning all of them and finding the geometry with he lowest RMSD with respect to the ensemble. 

Results are ambiguous as the average is not sufficient to evaluate variability among side chains and backbone, but for longer residues different peack can be distinguished among the RMSD distribution with respect to the average. Such peaks correspond to sidechains in an elongated position vs more packed ones.

To run the program, load the .pdb file inside the dades folder and use dis.py. 

Finally, in the animations folder one can see how all aligned structures compare to a given average residue, which can be created by vmd by the vmd_comandos.tcl script.

