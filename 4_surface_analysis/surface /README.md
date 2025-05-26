Aquest programa calcula la superfície d'una proteïna submergida dins d'una membrana lipídaca.
Demana el nom d'una proteïna "p"

Requereix un fitxer p.pdb amb la superfície efectiva de cada àtom de la proteïna guardada al b_factor. Això es pot dur a terme amb Pymol a través del pdb fent 

load_p
set dot_solvent = 1 #Per calcular la SASA 
set dot_density = 4 #Controla la precissió amb la que es calcula la SASA amb el 
set solvent_radius = 1.4 #Fixa el radi del disolvent
h_add
cmd.get_area(load_b=1)
rebuild
save .pdb

Utilitza el fitxer .pdb amb les dades de superfície efectiva per àtom calculades amb pymol i guardades al b_factor, per altra banda, necessita un fitxer p.txt on s'indica per ordre els intervals de residuus de la proteïna que es troben insertats teòricament en la membrana.

Per altra banda, al fitxer config també es pot indicar si es vol realitzar el càlcul considerant únicament certs aminoàcids específics o només certs àtoms específics de cada residuu. Els càlculs de la SASA accessible a cada àtom tenen els paràmetres indicats. 
