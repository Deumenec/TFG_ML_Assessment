Using these conformation lists, the conformations predicted by AlphaFold are compared to see if there are "crazy" conformations really different from those exhibited by the amino acids.

The goal of this section is to improve the results obtained for "average" amino acids. To do this, protein conformations from the PDB are analyzed by setting a certain RMSD threshold. Based on this RMSD, a list of possible conformations is created, and for each of these conformations, they are compared with the AlphaFold-predicted conformations to check if there are "crazy" conformations that differ greatly from those experimentally determined.

All pdb structures to create the database are stored in the dades folder while analized and compared structures are stored 