This folder contains predictions for proteins with salts present were weird behaiviours by AF3 were observed.
fold_2025_05_17_salty_prediction contains the 4 peptide chains corresponding with a potassium chanel, together with three RNA segments, two chloride ions and two potassium chanels. We observe that the model is wrongly understanding metal placement inside of the chanel, and is placing two potassium ions extremly close one to the other.

This can be explained by the fact that in the PDB model for this structure, where a potassium ion is modeled twice. When submiting a prediction with even more than 4 potassium ions, all S1, S2, S3, S4 positions are occupied so no problem at all, it is a possible transition state even though it is not really stable


