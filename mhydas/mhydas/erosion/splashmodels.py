import numpy as np

from mhydas.mhydas.utilities import variablesdefinition

def mean_weight_diameter(precipitation, inflow_unit, Q_CALC_UNIT, unit_area, SDR_ajust_Q,
                          global_param_as_dict,local_param_as_dict,KE_nu, KE_couv, kteste):
    # function [SPLASH_CALC_UNIT_LISEM, CONC_SPLASH, Splash_Unit_LISEM, Splash_direct, Splash_indirect] = f_MHYDAS_UH_Splash_MWD_Erosion(precipitation, inflow_unit, Q_CALC_UNIT, unit_width, unit_length, SDR_ajust_Q, PARAM,KE_nu, KE_couv, kteste)
    #
    # % Calcul de la production de sédiments par splash arrivant aux tronçons de rigole d'UN MOTIF à chaque pas de temps avec les équations de LISEM (g/pas de temps)
    # % Auteurs: Silvio Gumiere
    # % Version : 15/01/2007
    # % Fichier: f_MHYDAS_UH_Splash_Erosion.m
    # %
    # % Cette fonction définit la quantité de sédiments détachée par effet splash (Kg/pas de temps)
    # % au sein des unités élémentaires d'UN SEUL MOTIF et effectivement transférée dans le réseau de rigoles
    # % on distingue la zone "RILL" et la zone "INTERRILL"
    # %
    # % Les notations & les détails des équations sont donnés dans :
    # % 1- De Roo A.P.J & Wesseling C.G, 1996. LISEM : A single-event physically based hydrological
    # %    and soil erosion model for drainage basins. I : Theory, input & output
    # % 2- Morgan et al, 1992. EUROSEM documentation manual.
    #
    # %*********************    SORTIES DE LA FONCTION    ***********************
    # % SPLASH_CALC_UNIT_LISEM : quantité de sédiments issue du splash  arrivant à chaque noeud (en kg/pas de temps)
    # % CONC_SPLASH            : concentration en sédiments au sein des tronçons (en kg/m3 ou g/l)
    #
    # %********************   ARGUMENTS DE LA FONCTION    ***********************
    # % precipitation.mmh     : intensité de pluie au pas de temps 'dt' choisi (mm/h)
    # % precipitation.int     : quantité de pluie tombée (m/pas de temps)
    # % inflow_unit : débit d'entrée latéral (m3/s)
    # % Q_CALC_UNIT   : Débit d'entrée dans chaque noeud (m3/s); comprend entrée latérale et entrée par tronçon amont
    # % unit_width     : largeur d'une unité élémentaire (m)
    # % unit_length     : longueur d'une unité élémentaire (m)
    # % SDR_ajust_Q   : SDR affecté à chaque unité élémentaire en fonction du débit
    # % MWD           : Mean Weigth Diameter em (m) adaptée pour Silvio
    # % dt            : pas de temps (s) --> PARAM.dt
    # % Surf_couvert  : Surface du motif en % sous couvert végétal (% de la surface totale du motif) --> PARAM.Larg_couvert
    # % Haut_canopee  : Hauteur du couvert végétal occupant un motif --> PARAM.Haut_canopee
    # % KE            : Knetic Energy of Rainfall


    SPLASH_CALC_UNIT_LISEM = []
    CONC_SPLASH = []
    precipitation_values = precipitation[variablesdefinition.precipitation_label_custom]
    inflow_values = inflow_unit[variablesdefinition.streamflow_label_custom].values
    for i in range(int(local_param_as_dict[variablesdefinition.nb_unit])):# % garde fou : verifier Nb_unit > 1
        Splash_Unit_LISEM = []
        Splash_direct = []
        Splash_indirect = []
        intermed_Conc_Splash = []
        gap_counter = 0
        for j in range(len(precipitation)): #looks like the same is called for each elementary unit
            if precipitation_values[j] == 0 or Q_CALC_UNIT[j][i] == 0: #% Conditions: pluie et débit sur chaque unité élémentaire non nuls
                Splash_Unit_LISEM.insert(j, 0)
            else:
                #% Calcul du splash en Kg/pas de temps
                if j - gap_counter > 1:
                    Splash_direct.extend([0] * (j - gap_counter - 1))
                    Splash_indirect.extend([0] * (j - gap_counter - 1))

                _splash_direct = (2.1e-4/local_param_as_dict[variablesdefinition.mwd] * kteste) * KE_nu[j] * \
                                 precipitation_values[j] * \
                                        (1-local_param_as_dict[variablesdefinition.surf_couvert]) * \
                                        unit_area

                Splash_direct.append(_splash_direct)
                _splash_indirect = (2.1e-4/local_param_as_dict[variablesdefinition.mwd])*KE_couv * \
                                       precipitation_values[j] * \
                        local_param_as_dict[variablesdefinition.surf_couvert]*unit_area

                Splash_indirect.append(_splash_indirect)

                Splash_Unit_LISEM.append((_splash_direct + _splash_indirect) * SDR_ajust_Q[j])
                gap_counter = j

                #% Calcul des concentrations issues du splash
            _inflow = inflow_values[j]
            if _inflow == 0:
                intermed_Conc_Splash.append(0)
            else:
                intermed_Conc_Splash.append(Splash_Unit_LISEM[j]/
                                            (_inflow * global_param_as_dict[variablesdefinition.dt])
                               )
        #% Ecriture matrices
        if i == 0:
            SPLASH_CALC_UNIT_LISEM = np.c_[np.array(list(map(lambda x: x*0.5, Splash_Unit_LISEM)))]
            CONC_SPLASH = np.c_[np.array(intermed_Conc_Splash)]
        else:
            SPLASH_CALC_UNIT_LISEM = np.c_[SPLASH_CALC_UNIT_LISEM,
                                           np.array(list(map(lambda x: x*0.5, np.add(Splash_Unit_LISEM, MEM))))]
            CONC_SPLASH = np.c_[CONC_SPLASH, np.array(intermed_Conc_Splash)]

        MEM = Splash_Unit_LISEM

    SPLASH_CALC_UNIT_LISEM = np.c_[SPLASH_CALC_UNIT_LISEM, np.array(list(map(lambda x: x*0.5, Splash_Unit_LISEM)))]
    CONC_SPLASH = np.c_[CONC_SPLASH, np.array(intermed_Conc_Splash)]

    return SPLASH_CALC_UNIT_LISEM, CONC_SPLASH, Splash_Unit_LISEM, Splash_direct, Splash_indirect
