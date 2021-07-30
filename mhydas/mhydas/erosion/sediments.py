import math

from mhydas.mhydas.utilities import variablesdefinition


def adjusted_sediment_discharge_ratio(global_parameters_as_dict, local_parameters_as_dict, net_rainfall_as_array):
    # function [SDR_ajust_Q] = f_MHYDAS_UH_SDR_ajust_Q (net_precipitation, PARAM)
    #
    # % Calcul du SDR de chaque tronçon en fonction du débit issu du ruissellement sur la zone inter-rigole.
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Fichier: f_MHYDAS_UH_SDR_Pente_unit.m
    # %
    # % Le SDR est asservi au débit spécifique ruisselant (pluie nette) et est borné par une valeur maximale
    # % Les travaux réalisés sur l'érosion diffuse (Leguédois, 2003; Nord, 2006)
    # % donnent des résultats qui permettent un ajustement par un modèle exponentiel :
    # % SDR_ajust_Q = SDR_init(1-exp(-2.5*net_precipitation/Ruissel_max))
    # % D'après cette équation, quand net_precipitation/Ruissel_max = 1 ---> SDR = 0.91*SDR_max
    #
    # % *******************   ARGUMENTS DE LA FONCTION    ***********************
    # % net_precipitation          : pluie nette (mm / pas de temps)
    # % SDR_max     : Valeur du SDR pour un débit max (%) --> PARAM.SDR_max
    # % Ruissel_max :  Valeur du ruissellement nécessaire pour atteindre 91% du SDR_init (mm/pas de temps) --> PARAM.Ruissel_max

    adjusted_sediment_discharge_ratio = list(map(lambda x: local_parameters_as_dict[variablesdefinition.sdr_max] *
                                            (1-math.exp(-((2.5*3600000*x /
                                                         global_parameters_as_dict[variablesdefinition.dt]) /
                                            local_parameters_as_dict[variablesdefinition.ruissel_max]))),
                      net_rainfall_as_array))
    # ((PARAM.Larg_UH/PARAM.Nb_motifs)/local_parameters_as_dict[variablesdefinition.larg_rill]);
    # % Colocado um fator de distancia entre o rigole e a zona interrill mesmo que WEPP Silvio 20/03/2007
    return adjusted_sediment_discharge_ratio
    # %SDR_ajust_Q = PARAM.SDR_max*(1-exp(-(2.5*3600000/PARAM.Ruissel_max)*net_precipitation/global_parameters_as_dict[variablesdefinition.dt]))*
    # exp(-((PARAM.Larg_UH/PARAM.Nb_motifs)/local_parameters_as_dict[variablesdefinition.larg_rill]*0.000004))
    # % Colocado um fator de distancia entre o rigole e a zona interrill mesmo que WEPP Silvio 20/03/2007
    #  %   end
    #   %SDR_AJUST_Q(j,i) = SDR_ajust_Q;
    # %end


def sediment_production_per_time_interval(Masse_antecedent, Apport_Splash, Apport_Amont, Perte_Aval, conc_antecedent,
                                          pente_Tr, Vitesse_Eau_Tr, Long_unit, cc, dd, betha, global_parameters_as_dict, local_parameters_as_dict, kteste):
    #function [unit_flow_lisem, tc_lisem] = f_MHYDAS_UH_Flow_LISEM_Erosion(Masse_antecedent, Apport_Splash, Apport_Amont, Perte_Aval, conc_antecedent, pente_Tr, Vitesse_Eau_Tr, Long_unit, cc, dd, betha, PARAM, kteste)
    # % Calcul du dépôt OU de la production de sédiments par le ruissellement à chaque pas de temps
    # % et sur chaque unité élémentaire avec les équations de LISEM (kg/pas de temps)
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Fichier: f_MHYDAS_UH_Splash_Erosion.m
    #
    # % Cette fonction définit la quantité de sédiments détachée par le ruissellement (en Kg/pas de temps)
    # % lorsque la puissance du ruissellement le permet et lorsque la concentration en MES est inférieure
    # % à la capacité de transport.
    # % Si l'une, l'autre, ou les deux conditions ne sont pas respectées; alors la fonction définit la quantité de sédiment déposée
    #
    # % Les notations & les détails des équations sont     donnés dans :
    # % 1- De Roo A.P.J & Wesseling C.G, 1996. LISEM : A single-event physically based hydrological
    # %    and soil erosion model for drainage basins. I : Theory, input & output.
    # %    Hydrological Processes 10(8) :1107-1117.
    # % 2- Morgan R.P.C. et al, 1992. EUROSEM documentation manual, 34p.
    # % 3- Morgan R.P.C. et al, 1998. The European Soil Erosion Model (EUROSEM) : a dynamic approach for predicting
    # %    sediment transport from fields and small catchments. Earth Surf. Process. Landforms, 23 :527-544.
    # % 4- Soulsby, R.L., 1997. Dynamics of marine sands, a manual for practical applications. Thomas Telford
    # % 5- Guillaume Nord, 2006. (Thèse) Modélisation à bases phyiques des processus de l'érosion
    # %    à l'échelle de la parcelle --> VITESSE DE DEPOT (Soulsby)
    # % 6- Cadiergue et al, 1999. Vitesse de chute d'une particule lourde isolée dans un écoulement turbulent
    #
    # %    major changes :
    #
    # %    ver. 0.2: (16 octobre 2006) : Application de la loi de Soulsby (vitesse de dépôt )
    # %    --> prise en compte de la concentration dans le calcul de la vitesse de dépôt (introduction d'une dynamique temporelle)
    # %    --> utilisation du paramètre phi : degré d'homogénéité du profil vertical de concentration dan la lame d'eau
    # %    --> changement du paramètre viscosité cinématique (passage des poises aux m²/s)
    #
    # %*********************    SORTIES DE LA FONCTION    ***********************
    # % unit_flow_lisem  : Correspond à la production interne (kg/pas de temps)
    # % >0 si DETACHEMENT  / <0 si DEPOT
    #
    # %********************   ARGUMENTS DE LA FONCTION    ***********************
    # % Masse_antecedent : Masse de sédiments présents dans Tr au pas de tps précédent (kg)
    # % Apport_Splash    : Apport par « splash » en kg entre t-1 et t
    # % Apport_Amont     : masse de sédiment provenant du tronçon amont entre t-1 et t (kg)
    # % Perte_Aval       : masse de sédiment sortant du tronçon entre t-1 et t (kg)
    # % conc_antecedent  : Concentration en sédiments au pas de temps précedent (g/l ou kg/m3)
    # % pente_Tr         : pente (m/m) --> provient du vecteur pente_unit
    # % Hauteur_Eau_Tr   : Hauteur d'eau dans le tronçon au pas de temps considéré (m)
    # % Vitesse_Eau_Tr   : Vitesse de l'eau dans le tronçon au pas de temps considéré (m/s)
    # % Debit_Tr         : Débit d'entrée dans le tronçon au pas de temps considéré (m3/s)
    # % Long_unit        : longueur d'une unité élémentaire (m)
    # %-----------    PARAM   ---------------
    # % AGGRSTAB_LISEM   : Paramètre  d'érodibilité "stabilité structurale" --> nb de gouttes pour diminuer de moitié la taille des agrégats
    # % Dens_Sed         : Densité des sédiments  --> local_parameters_as_dict[variablesdefinition.dens_ed]
    # % D50_Sed          : Diamètre médian des particules (m) --> local_parameters_as_dict[variablesdefinition.d50_ed]
    # % Larg_rill        : largeur de la rigole supposée rectangulaire (m) --> local_parameters_as_dict[variablesdefinition.larg_rill]
    # % Dens_fluide      : Densité de l'eau  --> local_parameters_as_dict[variablesdefinition.dens_fluide]
    # % Visc_fluide      : Viscosité moléculaire cinématique de l'eau (m²/s) -->  local_parameters_as_dict[variablesdefinition.visc_fluide]
    # % Coeff_Phi        : Dégré d'homogénéité du profil vertical de concentration dans la lame d'eau [0-1]
    # % dt               : pas de temps (s) --> global_parameters_as_dict[variablesdefinition.dt]
    #
    # % **************** PROCESSUS DE DETACHEMENT ET DE DEPOT *******************
    # % ETAPE 1 : La capacité de transport est-elle positive?
    # %     i.e   : on compare le STREAM POWER à une valeur de STREAM POWER critique = à 0.004 m/s (Govers, 1990)
    # %     i.e   : Le ruissellement est-il suffisamment "puissant" pour transporter les particules?
    # %     Si non --> dépôt = w*C*Vs*dx*dt      si oui --> calcul de TC
    #
    # % ETAPE 2 : Si TC>0 --> Comparaison entre la capacité de transport (TC) et la charge en sédiments (C)
    # %    CAS 1  : TC<C --> Dépot (Tc-C)(1-exp(dt Vs/d))*h*w*dx
    # %    CAS 2  : TC>C --> Détachement
    #

    tc_lisem = 0#;
    #% 1) Comparaison entre Stream Power (w) & Critical Stream Power (wc)
    Delta_Stream_Pow = ((Vitesse_Eau_Tr *100 )*(pente_Tr))-0.4#; % De roo and Wesseling - HP - 1996(10)1107-1117

    #% Calcul du diamètre spécifique de la particule (Loi de Soulsby)
    D_specifique = ((((((local_parameters_as_dict[variablesdefinition.dens_ed]-local_parameters_as_dict[variablesdefinition.dens_fluide])/local_parameters_as_dict[variablesdefinition.dens_fluide])*9.81))/((local_parameters_as_dict[variablesdefinition.visc_fluide]))^2)^(1/3))*local_parameters_as_dict[variablesdefinition.d50_ed]#   ;
    #% Application de la formule de Soulsby
    Vit_depot = (local_parameters_as_dict[variablesdefinition.visc_fluide]/local_parameters_as_dict[variablesdefinition.d50_ed])*((10.36^2+1.049*(1-(conc_antecedent)/2650)^4.7*D_specifique^3)^(1/2)-10.36)#;
    #%Vit_depot = Vit_depot * PARAM.Coeff_Phi;

    if  Delta_Stream_Pow < 0:# % Puissance de ruissellement insuffisante pour transporter les particules
        unit_flow_lisem = -(local_parameters_as_dict[variablesdefinition.larg_rill]*conc_antecedent*Vit_depot*Long_unit).*global_parameters_as_dict[variablesdefinition.dt]#; % DEPOT : w*C*Vs*dx*dt (dépot en Kg/pas de temps)
        if unit_flow_lisem < -Masse_antecedent - Apport_Splash - Apport_Amont + Perte_Aval:
            unit_flow_lisem = -Masse_antecedent - Apport_Splash - Apport_Amont + Perte_Aval#;
    else:
        tc_lisem = local_parameters_as_dict[variablesdefinition.dens_ed]*1000*cc*math.pow(Delta_Stream_Pow,dd)#; % converti de m3/m3 en kg/m3 (prise en cpte de la masse volumique (2650kg/m3)
        #tc_lisem = TC#;
        #% 2) Comparaison entre capacité de transport et charge en sédiments
        if tc_lisem-conc_antecedent < 0 :#% Incapacité de transporter car concentration > capacité de transport --> DEPOT
            #% Equation de dépôt pour la version de LISEM postérieure à 2.02
            #% unit_flow_lisem = (TC-conc_antecedent)*(1-exp(global_parameters_as_dict[variablesdefinition.dt].*Vit_depot/Hauteur_Eau_Tr)).*Hauteur_Eau_Tr*local_parameters_as_dict[variablesdefinition.larg_rill]*Long_unit.*global_parameters_as_dict[variablesdefinition.dt]; % DEPOT : (Tc-C)(1-exp(dt Vs/d)*h*w*dx

            #% Equation de dépôt pour les version de LISEM 2.02 et antérieures
            unit_flow_lisem = (tc_lisem-conc_antecedent)*Vit_depot*local_parameters_as_dict[variablesdefinition.larg_rill]*Long_unit*global_parameters_as_dict[variablesdefinition.dt]#; % DEPOT : (Tc-C)Vs*w*dx
            if unit_flow_lisem < -Masse_antecedent - Apport_Splash - Apport_Amont + Perte_Aval:
                unit_flow_lisem = -Masse_antecedent - Apport_Splash - Apport_Amont + Perte_Aval#;

        else:
            unit_flow_lisem = betha *(tc_lisem-conc_antecedent)*Vit_depot *local_parameters_as_dict[variablesdefinition.larg_rill] *Long_unit *global_parameters_as_dict[variablesdefinition.dt]#; % ARRACHEMENT
    return unit_flow_lisem, tc_lisem

def sediments_concentration_per_unit(SPLASH_CALC_UNIT_LISEM, H_CALC_UNIT,  V_CALC_UNIT, Q_CALC_UNIT,
                                                       Q_CALC_SORTIE, unit_length, unit_width,
                                                       unit_slope, net_precipitation,
                                                        global_parameters_as_dict, local_parameters_as_dict):
    # function [CALC_CONC_TR_LISEM, CALC_PROD_INTERNE, MASSE_SED, CALC_TC_LISEM, CALC_VOL_TR_LISEM] = f_MHYDAS_UH_Conc_Tr_LISEM_bilan(SPLASH_CALC_UNIT_LISEM, H_CALC_UNIT, V_CALC_UNIT, Q_CALC_UNIT, Q_CALC_SORTIE, unit_length, unit_width, unit_slope, net_precipitation, PARAM, kteste)
    # % Calcul de la concentration en sédiments dans chaque tronçon (Kg/m3) & à chaque pas de temps avec les équations de LISEM (Kg/m3)
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Fichier: f_MHYDAS_UH_Conc_Tr_LISEM
    # %
    # % Les notations & les détails des équations sont donnés dans :
    # % 1- De Roo A.P.J & Wesseling C.G, 1996. LISEM : A single-event physically based hydrological
    # %    and soil erosion model for drainage basins. I : Theory, input & output.
    # %    Hydrological Processes 10(8) :1107-1117.
    # % 2- Morgan R.P.C. et al, 1992. EUROSEM documentation manual, 34p.
    # % 3- Morgan R.P.C. et al, 1998. The European Soil Erosion Model (EUROSEM) : a dynamic approach for predicting
    # %    sediment transport from fields and small catchments. Earth Surf. Process. Landforms, 23 :527-544.
    # % 4- Soulsby, R.L., 1997. Dynamics of marine sands, a manual for practical applications. Thomas Telford
    # % 5- Guillaume Nord, 2006. (Thèse) Modélisation à bases phyiques des processus de l'érosion
    # %    à l'échelle de la parcelle --> VITESSE DE DEPOT (Soulsby)
    # % 6- Cadiergue et al, 1999. Vitesse de chute d'une particule lourde isolée dans un écoulement turbulent
    #
    # %*********************    SORTIES DE LA FONCTION    ***********************
    # % CALC_CONC_TR_LISEM     : Concentration en sédiments au sein de chaque tronçon d'un motif & en SORTIE de motif (Kg/m3 <---> g/l)
    # % CALC_PROD_INTERNE      : Quantité de sédiments déposés (<0) ou érodés (>0) dans chaque tronçon (Kg/pas de temps)
    # % MASSE_SED              : Masse de sédiments en suspension au sein de chaque tronçon (Kg)
    # % CALC_TC_LISEM          : capacité de transport dans chaque tronçon en kg/m3 (g/l)
    # % CALC_VOL_TR_LISEM      : volume d'eau par tronçon (m3)
    #
    # %********************   ARGUMENTS DE LA FONCTION    ***********************
    # % SPLASH_CALC_UNIT_LISEM : quantité de sédiments issue du splash arrivant à chaque noeud (kg/pas de temps)
    # % H_CALC_UNIT            : Hauteur d'eau dans le tronçon (m)
    # % V_CALC_UNIT            : Vitesse de l'eau dans le tronçon (m/s)
    # % Q_CALC_UNIT            : Débit d'entrée dans le tronçon (m3/s)
    # % Q_CALC_SORTIE          : Débit de sortie dans le tronçon (m3/s)
    # % unit_length              : longueur d'une unité élémentaire (m)
    # % unit_width              : largeur d'une unité élémentaire (m)
    # % unit_slope             : pente des tronçons (m/m)
    # % net_precipitation                     : pluie nette (mm/pas de temps)
    # % COH                    : paramètre cohésion du sol de LISEM (déterminé à l'aide de tests Torvane) (kPa) --> local_parameters_as_dict[variablesdefinition.coh]
    # % Larg_rill              : largeur de la rigole supposée rectangulaire (m) --> local_parameters_as_dict[variablesdefinition.larg_rill]
    # % dt                     : pas de temps de calcul (s) --> global_parameters_as_dict[variablesdefinition.dt]
    #
    # %***********************   FONCTIONS APPELEES    **************************
    # % unit_flow_lisem        : DEPOT OU PRODUCTION de sédiments par le ruissellement à chaque pas de temps (Kg/pas de temps)
    #
    # %*************************   INITIALISATION    ****************************
    # %bilan volume
    CALC_VOL_TR_LISEM = []
    CALC_CONC_TR_LISEM = []
    CALC_PROD_INTERNE = []
    MASSE_SED = []
    CALC_TC_LISEM = []
    CALC_VOL_TR_LISEM.insert(0, [0]*local_parameters_as_dict[variablesdefinition.nb_unit])
    #%bilan masse sédiments
    CALC_CONC_TR_LISEM.insert(0, [0]*(local_parameters_as_dict[variablesdefinition.nb_unit] + 1))
    CALC_PROD_INTERNE.insert(0, [0]*local_parameters_as_dict[variablesdefinition.nb_unit])
    MASSE_SED.insert(0, [0]*local_parameters_as_dict[variablesdefinition.nb_unit])
    CALC_TC_LISEM.insert(0, [0]*local_parameters_as_dict[variablesdefinition.nb_unit])

    # % ******   CALCUL DE PARAMETRES INTERMEDIAIRES CONSTANTS POUR PRODUCTION INTERNE EN SEDIMENTS   *****
    # % fonction f_MHYDAS_UH_Flow_LISEM_Erosion
    # % Calcul des coefficients expérimentaux pour l'obtention de la capacité de transport
    cc = ((local_parameters_as_dict[variablesdefinition.d50_sed]*1000000 + 5) /0.32)**(-0.6)#; % ATTENTION : équation valide pour un diamètre de particules > 32 µm
    dd = ((local_parameters_as_dict[variablesdefinition.d50_sed]*1000000 + 5) /300)**0.25#;  % ATTENTION : équation valide pour un diamètre de particules > 32 µm
    # % Calcul du facteur d'efficacité du détachement lié à la cohésion du sol (Raws and Govers, 1988)
    # %YY = 1/(0.89+0.56.*local_parameters_as_dict[variablesdefinition.coh]);
    betha = 0.79*math.exp(-0.85*local_parameters_as_dict[variablesdefinition.coh])

    for j in range(1, len(Q_CALC_UNIT)):
        for i in range(local_parameters_as_dict[variablesdefinition.nb_unit]):#% garde fou : verifier Nb_unit > 1
            #%*********************     BILAN  VOLUME D'EAU (m3)   *********************
            if Q_CALC_UNIT[j][i] < 0.0000001 or H_CALC_UNIT[j][i] == 0:
                Volume_Tr_LISEM = 0
            else:
                #% Apport d'eau latéral par ruissellement sur zone diffuse en m3/pas de temps --> Q_entree_unit
                if i == 1:
                    Volume_Splash = 0.5*net_precipitation[j-1] * unit_length*unit_width
                    Volume_Amont = 0  # ;
                else:
                    Volume_Splash = net_precipitation[j-1] * unit_length*unit_width
                    Volume_Amont = Q_CALC_SORTIE[j-1][i-1] * global_parameters_as_dict[variablesdefinition.dt]
                #end
                 #% Apport d'eau à l'amont du tronçon (débit de sortie du tronçon précédent) en m3/pas de temps --> Q_CALC_SORTIE (t-1,Tr-1)
                # if i == 1:
                #    Volume_Amont = 0#;
                # else:
                #    Volume_Amont = Q_CALC_SORTIE(j-1,i-1)*global_parameters_as_dict[variablesdefinition.dt]#; %
                # #end
                 #% Sortie d'eau à l'aval du tronçon en m3/pas de temps --> Q_CALC_SORTIE (t-1,Tr)
                Volume_Aval = Q_CALC_SORTIE[j-1][i]*global_parameters_as_dict[variablesdefinition.dt]#;
                 #% volume d'eau dans le tronçon Tr au pas de temps précédent en m3
                Volume_antecedent = CALC_VOL_TR_LISEM[j-1][i]#;
                #% ECRITURE DU BILAN
                Volume_Tr_LISEM = Volume_antecedent  + Volume_Splash + Volume_Amont - Volume_Aval#;
                if Volume_Tr_LISEM < 0:
                    Volume_Tr_LISEM = 0#;

            #% Ecriture Matrice (stockage)
            CALC_VOL_TR_LISEM[j].insert(i, Volume_Tr_LISEM)#;

            #%**********************  BILAN DE MASSE  EN KG/dt   ***********************
            #% Le bilan de masse est réalisé en Kg/pas de temps
            if Q_CALC_UNIT[j][i] < 0.0000001 or H_CALC_UNIT[j][i] == 0:
                 Conc_Tr_LISEM = 0#;
                 Masse_Tr_LISEM = 0#;
                 CALC_TC_LISEM[j].insert(i, 0)#;
                 Prod_interne = 0#;
            else:
                #% Masse de sédiments dans le tronçon Tr au pas de temps précédent en Kg
                 Masse_antecedent = CALC_CONC_TR_LISEM[j-1][i]*CALC_VOL_TR_LISEM[j-1][i]#;
                #% Apport de sédiments par splash (au noeud amont) en Kg/pas de temps --> Q_entree_unit
                 Apport_Splash = SPLASH_CALC_UNIT_LISEM[j-1][i]#;
                #% Apport de sédiments à l'amont du tronçon (débit de sortie du tronçon précédent) en Kg/pas de temps --> Q_CALC_SORTIE (t-1,Tr-1)
                 if i == 1:
                     Apport_Amont = 0#;
                 else:
                     Apport_Amont = CALC_CONC_TR_LISEM[j-1][i-1]*Q_CALC_SORTIE[j-1][i-1]*global_parameters_as_dict[variablesdefinition.dt]#; %
                 #end
                #% Sortie de sédiments à l'aval du tronçon en Kg/pas de temps --> Q_CALC_SORTIE (t-1,Tr)
                 Perte_Aval = CALC_CONC_TR_LISEM[j-1][i]*Q_CALC_SORTIE[j-1][i]*global_parameters_as_dict[variablesdefinition.dt]#;
                #% Production interne au tronçon : Dépôt OU Arrachement par le ruissellement(en Kg/pas de temps) --> Q_CALC_UNIT
                 Prod_interne,tc_lisem = sediment_production_per_time_interval(Masse_antecedent, Apport_Splash, Apport_Amont, Perte_Aval, CALC_CONC_TR_LISEM[j-1][i], unit_slope[i], V_CALC_UNIT[j-1][i], unit_length, cc, dd, betha, PARAM)#;
                #%  Prod_interne = 0;
                #% TC = 0;
                 CALC_TC_LISEM[j][i] = tc_lisem#;
                #% ECRITURE DU BILAN
                 Masse_Tr_LISEM = Masse_antecedent + Prod_interne + Apport_Splash + Apport_Amont - Perte_Aval#;
                    #%if Masse_Tr_LISEM < 0 ;
                    #%   Masse_Tr_LISEM = 0;
                    #%end
                 if CALC_VOL_TR_LISEM[j][i] == 0:#;
                     Conc_Tr_LISEM = 0#;
                 else:
                     Conc_Tr_LISEM = Masse_Tr_LISEM / CALC_VOL_TR_LISEM[j][i]#;
                #end
            #end
            #% Ecriture Matrice (stockage)
            CALC_CONC_TR_LISEM[j].insert(i, Conc_Tr_LISEM)#;
            CALC_PROD_INTERNE[j].insert(i, Prod_interne)#;
            MASSE_SED[j].insert(i, Masse_Tr_LISEM)#;
            MASSE_SED[j].insert(i, Masse_Tr_LISEM)#;
        #end
        #% calcul de la concentration sortant du dernier noeud (colonne Nb_unit+1)
        if Q_CALC_UNIT[j][local_parameters_as_dict[variablesdefinition.nb_unit]+1] < 0.0000001 or Q_CALC_SORTIE[j][local_parameters_as_dict[variablesdefinition.nb_unit]] < 0.00000001:
            Conc_Tr_LISEM_aval = 0#;
            Masse_Tr_LISEM_aval = 0#;
        else:
            #% conservation de la matière (flux additif)
            Conc_Tr_LISEM_aval = (CALC_CONC_TR_LISEM[j][local_parameters_as_dict[variablesdefinition.nb_unit]]*Q_CALC_SORTIE[j][local_parameters_as_dict[variablesdefinition.nb_unit]]+(SPLASH_CALC_UNIT_LISEM[j-1][local_parameters_as_dict[variablesdefinition.nb_unit]+1]/global_parameters_as_dict[variablesdefinition.dt]))/Q_CALC_UNIT[j][local_parameters_as_dict[variablesdefinition.nb_unit]+1]#;
            Masse_Tr_LISEM_aval = Conc_Tr_LISEM_aval * Q_CALC_UNIT[j][local_parameters_as_dict[variablesdefinition.nb_unit]+1]*global_parameters_as_dict[variablesdefinition.dt]

            #% Ecriture Matrice (stockage du dernier noeud  )
            CALC_CONC_TR_LISEM[j].insert(local_parameters_as_dict[variablesdefinition.nb_unit]+1, Conc_Tr_LISEM_aval)#;
            MASSE_SED[j].insert(local_parameters_as_dict[variablesdefinition.nb_unit]+1, Masse_Tr_LISEM_aval)#;
    return CALC_CONC_TR_LISEM, CALC_PROD_INTERNE, MASSE_SED, CALC_TC_LISEM, CALC_VOL_TR_LISEM