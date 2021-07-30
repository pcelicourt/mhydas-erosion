import os
from datetime import datetime
from statistics import mean

import pandas as pd

from mhydas.mhydas.utilities import filesmanager, variablesdefinition, slopes
from mhydas.erosion.datetimetransformation import disaggregate_date_time_from_minute_to_seconds
from mhydas.mhydas.erosion import infiltration
from mhydas.mhydas.erosion import routing, rillproperties, kineticenergy, sediments, splashmodels


default_data_file_dir = "../data/"


class Model():
    def __init__(self, data_directory=None, main_config_file_name=None):
        if data_directory:
            self.data_directory = data_directory
        else:
            self.data_directory = default_data_file_dir

        self.default_main_config_file_name = "config.ini"
        self.default_global_param_config_file_name = "global_parameters.ini"
        self.default_specific_param_config_file_name = "local_parameters.ini"
        self.parameters = {}
        self.kteste = 1

        if main_config_file_name:
            _main_config_file_path = os.path.join(self.data_directory,
                                                  main_config_file_name
                                                 )
            return self.set_config_files_paths(_main_config_file_path)
        else:
            return self.set_config_files_paths()

    def set_config_files_paths(self, main_config_file_path=None):
        if not main_config_file_path:
            _path = os.path.join(self.data_directory, self.default_main_config_file_name)
            self.main_config_file_path = self.check_file_existence(_path)
            print(self.main_config_file_path)
        else:
            self.main_config_file_path = self.check_file_existence(main_config_file_path)

        if self.main_config_file_path:
            self.main_config_file_content = filesmanager.read_main_config_file(self.main_config_file_path)
            self.global_param_config_file_path = self.main_config_file_content.get(variablesdefinition.global_param)
            self.local_param_config_file_path = self.main_config_file_content.get(variablesdefinition.local_param)
        else:
            filesmanager.file_not_found_error("main configuration file")

    def set_parameters(self):
        if self.global_param_config_file_path and self.local_param_config_file_path:
            self.set_global_parameters(self.global_param_config_file_path)
            self.set_local_parameters(self.local_param_config_file_path)
        else:
            filesmanager.file_not_found_error("main configuration file")

    def set_global_parameters(self, file_path):
        self.parameters[variablesdefinition.global_param] =\
            filesmanager.read_model_parameters_config_file(self.global_param_config_file_path)

    def set_local_parameters(self, file_path):
        self.parameters[variablesdefinition.local_param] =\
            filesmanager.read_model_parameters_config_file(self.local_param_config_file_path)

    def check_file_existence(self, file_path):
        return (lambda _path: _path if os.path.exists(_path) else None)(file_path)

    def get_main_config_file(self):
        try:
            return self.main_config_file_path #return the content or the path?
        except:
            filesmanager.file_not_found_error("global configuration file")

    def get_global_parameters(self):
        try:
            return self.parameters[variablesdefinition.global_param]
        except:
            filesmanager.file_not_found_error("global configuration file")

    def get_local_parameters(self):
        try:
            return self.parameters[variablesdefinition.local_param]
        except:
            filesmanager.file_not_found_error("local configuration file")

    def get_precipitation_data(self):
        path = self.main_config_file_content.get(variablesdefinition.precipitation)
        data = disaggregate_date_time_from_minute_to_seconds(
                                        self.main_config_file_content.get(variablesdefinition.precipitation),
                                        "\t",
                                        self.parameters[variablesdefinition.global_param][variablesdefinition.dt]
        )
        return data

    def get_streamflow_data(self):
        columns = [variablesdefinition.datetime, variablesdefinition.streamflow_label,
                   variablesdefinition.streamflow_label_custom
                ]
        streamflow_data = pd.DataFrame(columns=columns)
        data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.streamflow),
                                        sep="\t", skiprows=2
                           )
        for index, row in data.iterrows():
            values = [datetime(*(row.values[:6])), row.values[6], row.values[6] * 0.001]
            streamflow_data = streamflow_data.append(dict(zip(columns, values)), ignore_index=True)
        return streamflow_data

    def get_sediment_concentration_data(self):
        columns = [variablesdefinition.datetime, variablesdefinition.concentration_label]
        mes_concentration_data = pd.DataFrame(columns=columns)
        data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.mes_concentration),
                                        sep="\t"
                           )
        for index, row in data.iterrows():
            print(*(row.values[:5]))
            values = [datetime(*list(map(int, row.values[:5]))), row.values[6]]
            mes_concentration_data = mes_concentration_data.append(dict(zip(columns, values)), ignore_index=True)
        return mes_concentration_data

    def get_infiltration(self):
        precipitation_data = self.get_precipitation_data()
        #print(self.parameters)
        _global_parameters = self.get_global_parameters()
        if _global_parameters[variablesdefinition.code_production] == 1:
            tetaet = (_global_parameters[variablesdefinition.tetai] - _global_parameters[variablesdefinition.tetar]) / \
                     (_global_parameters[variablesdefinition.tetas] - _global_parameters[variablesdefinition.tetar])
            storage_and_suction_factor = _global_parameters[variablesdefinition.hc] * (1 - 1 * (tetaet**6)) * \
                                         (_global_parameters[variablesdefinition.tetas] -
                                          _global_parameters[variablesdefinition.tetai]
                                          )

            self.time_steps, self.infiltration_capacity, self.net_precipitation = infiltration.morelseytouxmethod(precipitation_data,
                                                                       storage_and_suction_factor,
                                                                       _global_parameters
                                                                       )
        else:
            pass


    def get_inflow_unit(self):
        inflow_unit = self.net_precipitation * (unit_length * unit_width /_global_parameters[variablesdefinition.dt])#; % Débit d'entrée latéral (m3/s)
        return inflow_unit

    def get_inflow_from_net_precipitation(self):
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        # unit_length = (_local_parameters[variablesdefinition.dist_1] + _local_parameters[variablesdefinition.dist_2] +
        #              _local_parameters[variablesdefinition.dist_3]) / _local_parameters[variablesdefinition.nb_unit]#; % Longueur de l'unité élémentaire
        # unit_width = _local_parameters[variablesdefinition.larg_uh] / _local_parameters[variablesdefinition.nb_motifs]#; % Largeur de l'unité élémentaire (= largeur motif élémentaire)
        unit_slope = slopes.vectorize(_local_parameters)
        unit_celerity = _global_parameters[variablesdefinition.celerite] * unit_slope / mean(unit_slope)
        if not self.net_precipitation:
            self.get_infiltration()
        #compute inflow_unit as a separate method
        inflow_unit = self.get_inflow_unit()#self.net_precipitation * (unit_length * unit_width /_global_parameters[variablesdefinition.dt])#; % Débit d'entrée latéral (m3/s)
        outflow_unit = []
        if _global_parameters[variablesdefinition.code_transfert] == 1:
            print('2. MHYDAS_UH : Fonction de transfert par le modèle de l''onde diffusante')
        else:
            print('2. MHYDAS_UH : Fonction de transfert par le schéma de Crank-Nicholson')

        Q_CALC_UNIT = []
        #% Hauteur (m) d'eau en entrée du tronçon
        H_CALC_UNIT = []
        #% Vitesse (m/s) en entrée du tronçon
        V_CALC_UNIT = []
        #% Débit (m3/s) sortant du tronçon considéré (débit d'entrée, Q_entree, transféré)
        Q_CALC_SORTIE = []

        for i in range(_local_parameters[variablesdefinition.nb_unit]):#% garde fou : verifier Nb_unit > 1
            celerite_bief = unit_celerity[i]
            hayami_core = routing.hayamimodel(_global_parameters,unit_length,celerite_bief)#; % Calcul du Noyau d'Hayami
            #% Calcul du débit entrant dans un tronçon : entrée latérale (ruissellement interrill) + sortie du tronçon amont
            if i == 1:
               inflow =  0.5 * inflow_unit;
            else:
               inflow =  inflow_unit + outflow_unit#;%
            H_entree = []
            V_entree = []
            #% Calcul des hauteurs (m) et vitesses (m/s) dans la rigole en entrée d'unité élémentaire
            for j in range(len(inflow)):
                if inflow[j] ==0:
                    H_entree.insert(j, 0)
                    V_entree.insert(j, 0)
                else:
                    rill_width = _local_parameters[variablesdefinition.larg_rill]
                    _h_entree = rillproperties.water_depth(
                        inflow[j], _local_parameters[variablesdefinition.rugo_strickler],
                        rill_width,
                        unit_slope[i], _local_parameters[variablesdefinition.seuil_conv]
                    )
                    rill_section = rill_width * _h_entree
                    _v_entree = rillproperties.water_speed(inflow[j], rill_section)

                    H_entree.insert(j, _h_entree)
                    V_entree.insert(j, _v_entree)

            if _global_parameters[variablesdefinition.code_transfert] == 1:
                outflow_unit = routing.hayamitransfer(inflow, hayami_core, _global_parameters[variablesdefinition.dt])
            # % ********************* Crank Nicholson (numérique) ***********************
            elif _global_parameters[variablesdefinition.code_transfert] == 2:
                outflow_unit = routing.crancknicholsontransfer(inflow, celerite_bief,
                                                               _global_parameters[variablesdefinition.sigma],
                                                               unit_length, _global_parameters[variablesdefinition.dt],
                                                               len(self.time_steps),
                                                               _global_parameters[variablesdefinition.nb_dx],
                                                               _global_parameters[variablesdefinition.nb_dt],
                                                               _global_parameters[variablesdefinition.nb_pas_fictifs]
                                                               )
            else:
                raise
            #% ********** Calcul des variables HYDRO (communs aux 2 méthodes) **********
            #% Ecriture matrices
            #% Débit (m3/s) entrant dans un tronçon : entrée latérale (ruissellement interrill) + sortie du tronçon amont
            Q_CALC_UNIT.insert(i, inflow)
            #% Hauteur (m) d'eau en entrée du tronçon
            H_CALC_UNIT.insert(i, H_entree)
            #% Vitesse (m/s) en entrée du tronçon
            V_CALC_UNIT.insert(i, V_entree)
            #% Débit (m3/s) sortant du tronçon considéré (débit d'entrée, Q_entree, transféré)
            Q_CALC_SORTIE.insert(i, outflow_unit)

            #% Débit (m3/s) sortant d'un motif
            Q_sortie_motif = (outflow_unit + 0.5*inflow_unit);
            #% Débit (m3/s) sortant de la parcelle
            Q_sortie_parcelle = (outflow_unit + 0.5*inflow_unit)* _local_parameters[variablesdefinition.nb_motifs]
            #% Calcul des hauteurs et vitesses dans la rigole en sortie d'un motif
        for j in range(len(inflow)):
            if Q_sortie_motif[j] == 0:
                H_sortie_motif.insert(j, 0)
                V_sortie_motif.insert(j, 0)
            else:
                rill_width = _local_parameters[variablesdefinition.larg_rill]
                _water_depth = rillproperties.water_depth(Q_sortie_motif[j],
                                                            _local_parameters[variablesdefinition.rugo_strickler],
                                                            rill_width,
                                                            unit_slope[_local_parameters[variablesdefinition.nb_unit]],
                                                            _local_parameters[variablesdefinition.seuil_conv]
                                                        )
                H_sortie_motif.insert(j, _water_depth)
                _section = _water_depth * rill_width
                V_sortie_motif[j] = rillproperties.water_speed(Q_sortie_motif[j], _section)
        #% Ajout dans les  matrices (colonne Nb_unit+1) des valeurs de Q,H,V, pour le dernier noeud
        Q_CALC_UNIT.insert(_local_parameters[variablesdefinition.nb_unit]+1, Q_sortie_motif)
        H_CALC_UNIT.insert(_local_parameters[variablesdefinition.nb_unit]+1, H_sortie_motif)#; % Calcul avec pente du dernier tronçon
        V_CALC_UNIT.insert(_local_parameters[variablesdefinition.nb_unit]+1, V_sortie_motif)#; % Calcul avec pente du dernier tronçon
        return Q_CALC_UNIT, H_CALC_UNIT, V_CALC_UNIT, Q_CALC_SORTIE


    def get_sediment_production(self):
        #if erosion == 1 :
        # %***************************  EROSION PAR SPLASH  *************************
        # %Escolha dos metodos de calc. da energia cinetica
        # %case 1  = metodo de Zanchi e Torri
        # %case 2  = metodo de Hudson
        # %case 3  = metodo de Wischmeir e Smith
        # %case 4  = metodo de Kinnel
        # %case 5  = metodo de Kinnel - Miami
        # %case 6  = metodo de Marshall e Palmer
        # %case 7  = metodo de Carter - USA
        # %case 8  = metodo de McGregor and Mutchler - Mississippi
        # %case 9  = metodo de Rosewell - Brisbane
        # %case 10 = metodo de Rosewell - Gunnedah
        # %case 11 = metodo de Rosewell - Melbourne
        # %case 12 = metodo de Rosewell - Cowra
        # %case 13 = metodo de Rosewell - Zimbabwe
        # %case 14 = metodo de Rosewell - Miami

        kinetic_method = 1
        kteste = 1
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()

        unit_length = (_local_parameters[variablesdefinition.dist_1] + _local_parameters[variablesdefinition.dist_2] +
                     _local_parameters[variablesdefinition.dist_3]) / _local_parameters[variablesdefinition.nb_unit]#; % Longueur de l'unité élémentaire
        unit_width = _local_parameters[variablesdefinition.larg_uh] / _local_parameters[variablesdefinition.nb_motifs]#; % Largeur de l'unité élémentaire (= largeur motif élémentaire)
        unit_slope = slopes.vectorize(_local_parameters)

        _precipitation = self.get_precipitation_data()
        ke_nu = kineticenergy.rain_on_bare_soil(_precipitation,_global_parameters,kinetic_method)
        ke_veg = kineticenergy.rain_on_vegetation(_local_parameters)
        SDR_ajust_Q = sediments.adjusted_sediment_discharge_ratio(_global_parameters, _local_parameters, net_rainfall)
        Q_CALC_UNIT, H_CALC_UNIT, V_CALC_UNIT, Q_CALC_SORTIE = self.get_inflow_from_net_precipitation()

        type_splash = 2
        precipitation_data = self.get_precipitation_data()
        inflow_unit = self.get_inflow_unit()
        # if type_splash ==  1:
        #     method = 'LISEM'
        #     print('3. MHYDAS_UH_EROSION : LISEM + SDR  ')
        #     SPLASH_CALC_UNIT_LISEM, CONC_SPLASH,  Splash_Unit_LISEM, Splash_direct, Splash_indirect = f_MHYDAS_UH_Splash_LISEM_Erosion(precipitation_data, Q_entree_unit, Q_CALC_UNIT, Larg_unit, Long_unit, SDR_ajust_Q, PARAM,KE_nu, KE_couv, kteste)
        if type_splash == 2:
            method = 'MWD + SDR'
            print('3. MHYDAS_UH_EROSION : Calcul du SPLASH par MWD avec ajout SDR')
            SPLASH_CALC_UNIT_LISEM, CONC_SPLASH,  Splash_Unit_LISEM, Splash_direct, Splash_indirect = \
                splashmodels.mean_weight_diameter(precipitation_data, inflow_unit, Q_CALC_UNIT, unit_width,
                                                  unit_length, SDR_ajust_Q, _global_parameters, _local_parameters,
                                                  ke_nu, ke_veg, kteste)
        # elif type_splash == 3:
        #     method = 'EUROSEM'
        #     print('3. MHYDAS_UH_EROSION : Calcul du SPLASH par EUROSEM ')
        #     SPLASH_CALC_UNIT_LISEM, CONC_SPLASH,  Splash_Unit_LISEM, Splash_direct, Splash_indirect = f_MHYDAS_UH_Splash_EUROSEM_Erosion(Pluie, Q_entree_unit, Q_CALC_UNIT, Larg_unit, Long_unit, SDR_ajust_Q, PARAM, pn,KE_nu, KE_couv)
        #        #%splash_mesure (ki, length(Pluie.int)) = (Splash_Unit_LISEM);


            #%disp ('3bis. MHYDAS_UH_EROSION : Calcul du SPLASH par WEPP');
        #%[SPLASH_CALC_UNIT_WEPP] = f_MHYDAS_UH_Splash_WEPP_Erosion(Pluie, PARAM.Argile, PARAM.Sable, PARAM.Sablefin, Larg_unit, Long_unit, pn, pente_unit, SDR_pente,PARAM.dt, PARAM);

        #%******************    CALCUL CHARGE SOLIDE DANS TRONCONS   ***************

        print('4. MHYDAS_UH_EROSION : Calcul des concentrations en MES dans les tronçons')#;

        CALC_CONC_TR_LISEM, CALC_PROD_INTERNE, MASSE_SED, CALC_TC_LISEM, CALC_VOL_TR_LISEM = \
            sediments.sediments_concentration_per_unit(SPLASH_CALC_UNIT_LISEM, H_CALC_UNIT,  V_CALC_UNIT, Q_CALC_UNIT,
                                                       Q_CALC_SORTIE, unit_length, unit_width,
                                                       unit_slope, self.net_precipitation,
                                                       _global_parameters, _local_parameters
                                                       )#; % DAMIEN

       #% Ja esta usando a equaçao de Soulsby  para velocidade de deposito

        # %D_specifique = ((((((PARAM.Dens_Sed-PARAM.Dens_fluide)/PARAM.Dens_fluide)*9.81))/(PARAM.Visc_fluide)^2)^(1/3))*PARAM.D50_Sed
        # %Vit_depot = (PARAM.Visc_fluide/PARAM.D50_Sed)*((10.36^2+1.049*D_specifique^3)^(1/2)-10.36)
        # %kteste
        # %ki

if __name__ == '__main__':
    new_model = Model(data_directory=default_data_file_dir)
    new_model.set_parameters()
    #print(new_model.get_precipitation_data())
    #print(new_model.get_streamflow_data())
    #print(new_model.get_sediment_concentration_data())
    new_model.get_infiltration()
    # new_model.set_global_parameters_config_file()
    # new_model.set_specific_parameters_config_file()
    # print(filesmanager.read_main_config_file(new_model.main_config_file_path))
    # print(filesmanager.read_model_parameters_config_file(new_model.global_param_config_file_path))
    # print(filesmanager.read_model_parameters_config_file(new_model.specific_param_config_file_path))


    # def set_config_file(self, file_path, default_path=None):
    #     if file_path:
    #         return self.check_file_existence(file_path)
    #     else:
    #         if default_path:
    #             return self.check_file_existence(default_path)
    #         else:
    #             return None
    #
    # def set_main_config_file(self, file_path=None):
    #     if file_path:
    #         self.main_config_file_path = self.set_config_file(file_path)
    #     else:
    #         default_config_file_path = glob.glob(os.path.join(self.data_directory,
    #                                                           self.default_main_config_file
    #                                                           )
    #                                              )[-1]
    #         self.main_config_file_path = self.set_config_file(default_config_file_path)
    #
    # def set_global_parameters_config_file(self, file_path=None):
    #     if file_path:
    #         self.global_param_config_file_path = self.set_config_file(file_path)
    #     else:
    #         default_config_file_path = glob.glob(os.path.join(self.data_directory,
    #                                                           self.default_global_param_config_file)
    #                                              )[-1]
    #         self.global_param_config_file_path = self.set_config_file(default_config_file_path)
    #
    # def set_specific_parameters_config_file(self, file_path=None):
    #     if file_path:
    #         self.specific_param_config_file_path = self.set_config_file(file_path)
    #     else:
    #         default_config_file_path = glob.glob(os.path.join(self.data_directory,
    #                                                           self.default_specific_param_config_file)
    #                                              )[-1]
    #         self.specific_param_config_file_path = self.set_config_file(default_config_file_path)