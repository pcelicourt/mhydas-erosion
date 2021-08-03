import os
from datetime import datetime
from statistics import mean
import operator

import pandas as pd
import numpy as np

from mhydas.mhydas.utilities import filesmanager, variablesdefinition, slopes, graphics
from mhydas.mhydas.erosion import infiltration
from mhydas.mhydas.erosion import routing, precipitation, rillproperties, kineticenergy, \
    sediments, splashmodels, errors


default_data_file_dir = "../data/"

#https://base-n.de/matlab/code_beautifier.html
class Model:
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
        self.net_precipitation = []

        if main_config_file_name:
            _main_config_file_path = os.path.join(self.data_directory,
                                                  main_config_file_name
                                                 )
            self.set_config_files_paths(_main_config_file_path)
        else:
            self.set_config_files_paths()

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
            self.set_global_parameters()
            self.set_local_parameters()
        else:
            filesmanager.file_not_found_error("main configuration file")

    def set_global_parameters(self):
        self.parameters[variablesdefinition.global_param] =\
            filesmanager.read_model_parameters_config_file(self.global_param_config_file_path)

    def set_local_parameters(self):
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
        data = precipitation.disaggregate_date_time_from_minute_to_seconds(
                                        self.main_config_file_content.get(variablesdefinition.precipitation),
                                        "\t",
                                        self.parameters[variablesdefinition.global_param][variablesdefinition.dt]
        )
        return data

    def get_streamflow_data(self):
        columns = [variablesdefinition.timestamp, variablesdefinition.streamflow_label,
                   variablesdefinition.streamflow_label_custom
                ]
        streamflow_data = pd.DataFrame(columns=columns)
        data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.streamflow),
                                        sep="\t", skiprows=2
                           )
        for index, row in data.iterrows():
            values = [datetime.timestamp(datetime(*(row.values[:6]))), row.values[6], row.values[6] * 0.001]
            streamflow_data = streamflow_data.append(dict(zip(columns, values)), ignore_index=True)
        return streamflow_data

    def get_sediment_concentration_data(self):
        columns = [variablesdefinition.timestamp, variablesdefinition.concentration_label]
        mes_concentration_data = pd.DataFrame(columns=columns)
        data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.mes_concentration),
                                        sep="\t"
                           )
        for index, row in data.iterrows():
            values = [datetime.timestamp(datetime(*list(map(int, row.values[:5])))), row.values[6]]
            mes_concentration_data = mes_concentration_data.append(dict(zip(columns, values)), ignore_index=True)
        return mes_concentration_data

    def get_infiltration_data(self):
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
            #TO DO: check this portion for the return values
            self.time_steps, self.infiltration_capacity, self.net_precipitation = \
                infiltration.morelseytouxmethod(precipitation_data, storage_and_suction_factor, _global_parameters)
            return self.time_steps, self.infiltration_capacity, self.net_precipitation
        else:
            pass


    def get_inflow_unit(self):
        _local_parameters = self.get_local_parameters()
        _global_parameters = self.get_global_parameters()
        unit_length = _local_parameters[variablesdefinition.long_uh]
        unit_width = _local_parameters[variablesdefinition.larg_uh]
        if not len(self.net_precipitation):
            self.get_infiltration_data()
        inflow_unit = self.net_precipitation * (unit_length * unit_width /_global_parameters[variablesdefinition.dt])#; % Débit d'entrée latéral (m3/s)
        return inflow_unit

    def get_inflow_from_net_precipitation(self):
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        unit_length = _local_parameters[variablesdefinition.long_uh]
        unit_width = _local_parameters[variablesdefinition.larg_uh]
        # unit_length = (_local_parameters[variablesdefinition.dist_1] + _local_parameters[variablesdefinition.dist_2] +
        #              _local_parameters[variablesdefinition.dist_3]) / _local_parameters[variablesdefinition.nb_unit]#; % Longueur de l'unité élémentaire
        # unit_width = _local_parameters[variablesdefinition.larg_uh] / _local_parameters[variablesdefinition.nb_motifs]#; % Largeur de l'unité élémentaire (= largeur motif élémentaire)
        unit_slope = slopes.vectorize(_local_parameters)
        unit_celerity = _global_parameters[variablesdefinition.celerite] * unit_slope / mean(unit_slope)
        if not self.net_precipitation:
            self.get_infiltration_data()
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
               inflow =  0.5 * inflow_unit
            else:
               inflow =  inflow_unit + outflow_unit
            H_entree = []
            V_entree = []
            #% Calcul des hauteurs (m) et vitesses (m/s) dans la rigole en entrée d'unité élémentaire
            for j in range(len(inflow)):
                if inflow[j] == 0:
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
        self.Q_sortie_motif = (outflow_unit + 0.5*inflow_unit)
        #% Débit (m3/s) sortant de la parcelle
        self.Q_sortie_parcelle = (outflow_unit + 0.5*inflow_unit)* _local_parameters[variablesdefinition.nb_motifs]
        #% Calcul des hauteurs et vitesses dans la rigole en sortie d'un motif
        H_sortie_motif = []
        V_sortie_motif = []
        #Q_sortie_motif = []
        for j in range(len(inflow)):
            if self.Q_sortie_motif[j] == 0:
                H_sortie_motif.insert(j, 0)
                V_sortie_motif.insert(j, 0)
            else:
                rill_width = _local_parameters[variablesdefinition.larg_rill]
                _water_depth = rillproperties.water_depth(self.Q_sortie_motif[j],
                                                            _local_parameters[variablesdefinition.rugo_strickler],
                                                            rill_width,
                                                            unit_slope[_local_parameters[variablesdefinition.nb_unit]],
                                                            _local_parameters[variablesdefinition.seuil_conv]
                                                        )
                H_sortie_motif.insert(j, _water_depth)
                _section = _water_depth * rill_width
                V_sortie_motif.insert(j, rillproperties.water_speed(self.Q_sortie_motif[j], _section))
        #% Ajout dans les  matrices (colonne Nb_unit+1) des valeurs de Q,H,V, pour le dernier noeud
        Q_CALC_UNIT.insert(_local_parameters[variablesdefinition.nb_unit]+1, self.Q_sortie_motif)
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
        #kteste = 1
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()

        unit_length = (_local_parameters[variablesdefinition.dist_1] + _local_parameters[variablesdefinition.dist_2] +
                     _local_parameters[variablesdefinition.dist_3]) / _local_parameters[variablesdefinition.nb_unit]#; % Longueur de l'unité élémentaire
        unit_width = _local_parameters[variablesdefinition.larg_uh] / _local_parameters[variablesdefinition.nb_motifs]#; % Largeur de l'unité élémentaire (= largeur motif élémentaire)
        unit_slope = slopes.vectorize(_local_parameters)

        precipitation_data = self.get_precipitation_data()
        #_precipitation = self.get_precipitation_data()
        ke_nu = kineticenergy.rain_on_bare_soil(precipitation_data,_global_parameters,kinetic_method)
        ke_veg = kineticenergy.rain_on_vegetation(_local_parameters)
        SDR_ajust_Q = sediments.adjusted_sediment_discharge_ratio(_global_parameters, _local_parameters, self.net_precipitation)
        self.Q_CALC_UNIT, H_CALC_UNIT, V_CALC_UNIT, Q_CALC_SORTIE = self.get_inflow_from_net_precipitation()

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
            self.SPLASH_CALC_UNIT_LISEM, CONC_SPLASH,  Splash_Unit_LISEM, self.Splash_direct, self.Splash_indirect = \
                splashmodels.mean_weight_diameter(precipitation_data, inflow_unit, self.Q_CALC_UNIT, unit_width,
                                                  unit_length, SDR_ajust_Q, _global_parameters, _local_parameters,
                                                  ke_nu, ke_veg, self.kteste)
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
            sediments.sediments_concentration_per_unit(self.SPLASH_CALC_UNIT_LISEM, H_CALC_UNIT,  V_CALC_UNIT,
                                                       self.Q_CALC_UNIT, Q_CALC_SORTIE, unit_length, unit_width,
                                                       unit_slope, self.net_precipitation,
                                                       _global_parameters, _local_parameters
                                                       )
        return CALC_CONC_TR_LISEM, CALC_PROD_INTERNE, MASSE_SED, CALC_TC_LISEM, CALC_VOL_TR_LISEM

       #% Ja esta usando a equaçao de Soulsby  para velocidade de deposito

        # %D_specifique = ((((((PARAM.Dens_Sed-PARAM.Dens_fluide)/PARAM.Dens_fluide)*9.81))/(PARAM.Visc_fluide)^2)^(1/3))*PARAM.D50_Sed
        # %Vit_depot = (PARAM.Visc_fluide/PARAM.D50_Sed)*((10.36^2+1.049*D_specifique^3)^(1/2)-10.36)
        # %kteste
        # %ki

    def get_hydrologic_balance(self):
        #%*************************   BILAN HYDROLOGIQUE   *************************
        #% Calcul des volumes : Pluie totale, Lame infiltrée, lame ruisselée, Volume mesuré, volume calculé
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        precipitation_data = self.get_precipitation_data()
        if not len(self.net_precipitation):
            self.get_infiltration_data()
        L_Pluie = sum(precipitation_data[variablesdefinition.precipitation_label_custom].values) * 1000
        L_Ruiss = sum(self.net_precipitation) * 1000
        L_Inf = L_Pluie - L_Ruiss#; % Lames précipitées, ruisselées et infiltrées en mm
        streamflow = self.get_streamflow_data()
        Vol_mes = np.trapz(streamflow[variablesdefinition.timestamp].values,
                           streamflow[variablesdefinition.streamflow_label_custom].values*86400)#; % Calcul du volume écoulé mesuré en m3


        Vol_cal = sum(self.Q_sortie_parcelle)*_global_parameters[variablesdefinition.dt]#; % Calcul du volume écoulé caclulé en m3
        Qmax_mes,imax = max(enumerate(streamflow[variablesdefinition.streamflow_label_custom].values),
                            key=operator.itemgetter(1))
        Qmax_cal,imax = max(enumerate(self.Q_sortie_parcelle), key=operator.itemgetter(1))

        aux = set(precipitation_data[variablesdefinition.timestamp].values).intersection(
                              streamflow[variablesdefinition.timestamp].values)
        ii = [precipitation_data[variablesdefinition.timestamp].values.index(x) for x in aux]
        jj = [streamflow[variablesdefinition.timestamp].values.index(x) for x in aux]
        coeff_Nash, RMSE, CRM = errors.nash_criteria(streamflow[variablesdefinition.streamflow_label_custom].values[jj],
                                                 self.Q_sortie_parcelle[ii])


    def get_erosion_balance(self):
        #%*******************   BILAN EROSION (mesurée & simulée) *******************
        #% Bilan au tronçon & à la parcelle
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        streamflow = self.get_streamflow_data()
        mes = self.get_sediment_concentration_data()
        CALC_CONC_TR_LISEM, CALC_PROD_INTERNE, MASSE_SED, CALC_TC_LISEM, CALC_VOL_TR_LISEM = \
        self.get_sediment_production()
        CALC_Prod_interne_Tr, CALC_Sortie_MES_Parcelle, CALC_Splash_Direct_Tot_Parcelle, \
        CALC_Splash_Indirect_Tot_Parcelle, CALC_Splash_Effectif_Parcelle = \
        sediments.sediments_balance(self.Q_CALC_UNIT,
                                  CALC_CONC_TR_LISEM,
                                  CALC_PROD_INTERNE,
                                  self.Splash_direct,
                                  self.Splash_indirect,
                                  self.SPLASH_CALC_UNIT_LISEM,
                                  _local_parameters,
                                  _global_parameters)
        #% Calcul des concentrations max (mesurée et simulée)
        Cmax_mes = max(mes[variablesdefinition.mes_concentration].values)
        Cmax_cal = max(CALC_CONC_TR_LISEM[:, _local_parameters[variablesdefinition.nb_unit]+1])
        #%calcul de la masse en sédiment mesurée exportée
        sed_mes = sediments.measured_sediment_mass(streamflow, mes)


    def create_sedimentograph(self):
        precipitation_data = self.get_precipitation_data()
        streamflow_data = self.get_streamflow_data()
        _, infiltration_data, _ = self.get_infiltration_data()
        mes = self.get_sediment_concentration_data()
        #print(precipitation_data.head())
        graphics.sedimentograph(precipitation_data, infiltration_data, streamflow_data, mes)



if __name__ == '__main__':
    new_model = Model(data_directory=default_data_file_dir)
    new_model.set_parameters()
    new_model.create_sedimentograph()
    #print(new_model.get_precipitation_data())
    #print(new_model.get_streamflow_data())
    #print(new_model.get_sediment_concentration_data())
    #new_model.get_infiltration_data()
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