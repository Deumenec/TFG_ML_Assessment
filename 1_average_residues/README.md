3/2/2025
Aquest és el primer programa del TFG, l'objectiu é validar que l'estructura dels aminoacids proporcionada per alphafold es correspon amb l'estructura que s'observa per a aquests als aminoacids normals.

Amb aquest objectiu, he escrit un programa que calcula la forma mitjana dels residuus corresponents a un aminoacid a una certa proteïna. El criteri per càlcular aquesta forma mitjana és determinar aquella geometria amb una mínima RMSD respecte la resta tots els residuus alineats amb aquesta. 

Els resultats obtinguts són ambiguus. En general però, s'observen diferències significatives entre les distribucions obtingudes. En general, els resultats d'alphafold semblen més centrats al voltant de certes conformacions per als aminoàcids mentres que d'altres no són tant presents.

Hipòtesis: Com alphafold està entrenat amb molècules que cristalitzen, té una certa tendència a replicar exactament aquestes estructures per a cada residuu en comptes de donar lloc a la fluxionalitat natural d'aquests. De tota manera, els resultats obtinguts.

De tota manera, els resultats obtinguts són massa generals. Amb l'objectiu de refinar-los, es vol fer un programa on l'estructura respecte la cual es comparen les estructures sigui una mica més fina.
