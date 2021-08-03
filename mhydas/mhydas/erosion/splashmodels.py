from mhydas.mhydas.utilities import variablesdefinition

def mean_weight_diameter(precipitation, inflow_unit, Q_CALC_UNIT, unit_width, unit_length, SDR_ajust_Q, local_param_as_dict, global_param_as_dict,KE_nu, KE_couv, kteste):
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
    Splash_Unit_LISEM = []
    Splash_direct = []
    Splash_indirect = []
    SPLASH_CALC_UNIT_LISEM = []
    CONC_SPLASH = []
    intermed_Conc_Splash = []
    for i in range(local_param_as_dict[variablesdefinition.nb_unit]):# % garde fou : verifier Nb_unit > 1
        for j in range(len(precipitation)):
            if precipitation[variablesdefinition.precipitation_label][j] == 0 or Q_CALC_UNIT[j][i] == 0: #% Conditions: pluie et débit sur chaque unité élémentaire non nuls
                Splash_Unit_LISEM.insert(j, 0)
            else:
                #% Calcul du splash en Kg/pas de temps
                Splash_direct.insert(j, (((2.1e-4/local_param_as_dict[variablesdefinition.mwd]*kteste)*KE_nu[j]*
                                      (precipitation[variablesdefinition.precipitation_label_custom][j]))*(
                                        (1-local_param_as_dict[
                                            variablesdefinition.surf_couvert
                                        ]
                                         )*unit_width*unit_length)
                                     )
                                  )
                Splash_indirect.insert(j, ((((2.1e-4/local_param_as_dict[variablesdefinition.mwd])*KE_couv))*
                                       (precipitation[variablesdefinition.precipitation_label_custom][j]))*(
                        local_param_as_dict[variablesdefinition.surf_couvert]*unit_width*unit_length))

                Splash_Unit_LISEM.insert(j, (Splash_direct[j] + Splash_indirect[j])*SDR_ajust_Q[j])

                #% Calcul des concentrations issues du splash
            if inflow_unit[j] == 0:
                intermed_Conc_Splash.insert(j, 0)
            else:
                intermed_Conc_Splash.insert(j, Splash_Unit_LISEM[j]/
                               (inflow_unit[j]*global_param_as_dict[variablesdefinition.dt])
                               )
       #% Ecriture matrices
        if i == 1:
            SPLASH_CALC_UNIT_LISEM.insert(i, 0.5*Splash_Unit_LISEM)
        else:
            SPLASH_CALC_UNIT_LISEM.insert(i, 0.5*Splash_Unit_LISEM+0.5*MEM)

        MEM = Splash_Unit_LISEM
        CONC_SPLASH.insert(i, intermed_Conc_Splash)

    SPLASH_CALC_UNIT_LISEM.insert(local_param_as_dict[variablesdefinition.nb_unit] + 1, 0.5*Splash_Unit_LISEM)
    CONC_SPLASH.insert(local_param_as_dict[variablesdefinition.nb_unit]+1, intermed_Conc_Splash)
    return SPLASH_CALC_UNIT_LISEM, CONC_SPLASH, Splash_Unit_LISEM, Splash_direct, Splash_indirect
