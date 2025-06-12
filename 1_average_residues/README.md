3/2/2025
This is the first program wrote in the TFG with the goal of validating aminoacid structure proportioned by AF3. It calculates the average shape for each aminoacid by aligning all of them and finding the geometry with he lowest RMSD with respect to the ensemble. 

Results are ambiguous as the average is not sufficient to evaluate variability among side chains and backbone, but for longer residues different peaks can be distinguished among the RMSD distribution with respect to the average. Such peaks correspond to sidechains in an elongated position vs more packed ones.

To run the program, load the .pdb file inside the dades folder and use dis.py. Aligned sequences of amino acids and average shapes for the given protein will be saved in their corresponding folder inside results.

"vmd_comandos.tcl" creates an animation using VMD of all found conformations when compared to the "average conformation" calculated using the calculated sequences.

11/6/2025
This section was not included in the full work owing to its similarity to the residue conformational analysis, which provided more significant results. 