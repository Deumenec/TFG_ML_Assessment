# VMD for MACOSXARM64, version 1.9.4a57 (April 27, 2022)
# Log file '/Users/deumenec/Documents/Uni/TFG química general/TFG química/TFG_ML_Assesment/distancia_atomica/dades/vmd_comandos', created by user deumenec
display resetview
display resetview
mol addrep 0
display resetview
mol new {/Users/deumenec/Documents/Uni/TFG química general/TFG química/TFG_ML_Assesment/distancia_atomica/dades/9kmh.cif} type {pdbx} first 0 last -1 step 1 waitfor 1
animate style Loop
menu save off
menu save on

animate write pdb {/Users/deumenec/Documents/Uni/TFG química general/TFG química/TFG_ML_Assesment/distancia_atomica/dades/9kmh.pdb} be g 0 end 0 skip 1 0
quit
