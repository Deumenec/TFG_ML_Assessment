L'entrada 4qry del PDB es tracta d'una pharaonis halorhodopsin en forma de trímer i complexada amb un io de brom. És tracta de tres alpha helixes insertades en una membrana que al seu torn poden formar un trímer tal i com es modela a l'estructura pujada. 
En general, la proteïna s'encarrega del bombeig d'halurs a través d'una membrana desencadenat per la presència de llum. 

Es genera una predicció d'Alphafold a través de la seqüència FASTA subministrada i després es realitzen les següents seleccions per considerar les mateixes parts de les molecules:

ORIGINAL: protein and chain A and altloc A and not residue 18

PREDITA: protein and residue 17 to 276 and not residue 18

L'estructura modelada per Alphafold presenta un TM score de 0.9905, igualment però, en pujar les dues estructures al PPM server, la predita per AF3 no presenta cap energia d'estabilització per trobar-se insertada en una membrana mentres que la original presenta una $\Delta$G de -73.6 kcal/mol.

A continuació, busquem una explicació raonable per a aquesta diferència substancial. A nivell de Ramachandran plot, no s'observa cap diferència significativa entre les dues estructures predites, només potser un major coeficient de correlació entre Psi i Phi a l'estructura experimental que al model, encara que segurament no prou significatiu.

Per altra banda, amb MDAnalysis, es prediu la SASA i la polar SASA per comprovar com canvien entre les dues estructures predites: 
