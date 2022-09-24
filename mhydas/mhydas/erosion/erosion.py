#import os
#from datetime import datetime
#from statistics import mean
import operator

import pandas as pd
import numpy as np

from mhydas.mhydas.utilities import variablesdefinition, graphics
from mhydas.mhydas.erosion import routing, rillproperties, kineticenergy, \
    sediments, splashmodels, errors, infiltration

default_data_file_dir = "../data/"


class Model:
    def __init__(self, global_parameters, local_parameters=None):
        self.local_parameters = local_parameters
        self.global_parameters = global_parameters
        self.parameters = {}
        self.kteste = 1
        self.net_precipitation = pd.DataFrame()

    def set_parameters(self):
        self.set_global_parameters()
        self.set_local_parameters()

    def set_global_parameters(self):
        self.parameters[variablesdefinition.global_param] = self.global_parameters

    def set_local_parameters(self):
        self.parameters[variablesdefinition.local_param] = self.local_parameters
            
    def get_global_parameters(self):
        try:
            return self.parameters[variablesdefinition.global_param]
        except KeyError:
            raise KeyError("Global parameters not set.")

    def get_local_parameters(self):
        try:
            return self.parameters[variablesdefinition.local_param]
        except KeyError:
            raise KeyError("Local parameters not set.")

    # def get_streamflow_data(self):
    #     columns = [variablesdefinition.datetime, variablesdefinition.streamflow_label,
    #                variablesdefinition.streamflow_label_custom
    #                ]
    #     streamflow_data = pd.DataFrame(columns=columns)
    #     data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.streamflow),
    #                        sep="\t", skiprows=2, header=None
    #                        ).convert_dtypes()
    #     for index, row in data.iterrows():
    #         _time = datetime(*(row.values[:6]))
    #         values = [_time, row.values[6], row.values[6] * 0.001]
    #         streamflow_data = streamflow_data.append(dict(zip(columns, values)), ignore_index=True)
    #     return streamflow_data

# =============================================================================
#     def get_sediment_concentration_data(self):
#         columns = [variablesdefinition.datetime, variablesdefinition.concentration_label] #variablesdefinition.timestamp,
#         mes_concentration_data = pd.DataFrame(columns=columns)
#         data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.mes_concentration),
#                            sep="\t", skiprows=3, header=None
#                            ).convert_dtypes()
#         for index, row in data.iterrows():
#             _time = datetime(*list(map(int, row.values[:6])))
#             values = [_time, row.values[6]] #datetime.toordinal(_time) + 366,
#             mes_concentration_data = mes_concentration_data.append(dict(zip(columns, values)), ignore_index=True)
#         #print("mes", mes_concentration_data.head())
#         return mes_concentration_data
# =============================================================================

    def get_infiltration_data(self):
        _global_parameters = self.get_global_parameters()
        if _global_parameters[variablesdefinition.code_production] == 1:
            tetaet = (_global_parameters[variablesdefinition.tetai] - _global_parameters[variablesdefinition.tetar]) / \
                     (_global_parameters[variablesdefinition.tetas] - _global_parameters[variablesdefinition.tetar])
            storage_and_suction_factor = _global_parameters[variablesdefinition.hc] * (1 - tetaet**6) * \
                                         (_global_parameters[variablesdefinition.tetas] -
                                          _global_parameters[variablesdefinition.tetai]
                                          )
            # TO DO: check this portion for the return values
            self.time_steps, self.infiltration_capacity, self.net_precipitation, self.infiltration_data = \
                infiltration.morelseytouxmethod(self.precipitation_data, storage_and_suction_factor, _global_parameters)
            return self.time_steps, self.infiltration_capacity, self.net_precipitation, self.infiltration_data
        else:
            pass

    def get_inflow_unit(self):
        _local_parameters = self.get_local_parameters()
        _global_parameters = self.get_global_parameters()
        unit_length = _local_parameters[variablesdefinition.long_uh]  
        unit_width = _local_parameters[variablesdefinition.larg_uh] / _local_parameters[variablesdefinition.nb_motifs]
        if not len(self.net_precipitation):
            self.get_infiltration_data()

        inflow_unit = pd.DataFrame({variablesdefinition.datetime: self.net_precipitation[variablesdefinition.datetime],
                                    variablesdefinition.streamflow_label_custom: list(map(lambda value: value *
                                                                                                       unit_length *
                                                                                                       unit_width /
                                                                        _global_parameters[variablesdefinition.dt],
                                            self.net_precipitation[variablesdefinition.precipitation_label_custom]))
                                    }
                                   )
        return inflow_unit

    def get_inflow_from_net_precipitation(self):
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        unit_length = _local_parameters[variablesdefinition.long_uh]
        # unit_length = (_local_parameters[variablesdefinition.dist_1] + _local_parameters[variablesdefinition.dist_2] +
        #                _local_parameters[variablesdefinition.dist_3]) / _local_parameters[variablesdefinition.nb_unit]
        unit_slope = _local_parameters[variablesdefinition.pente_uh]#slopes.vectorize(_local_parameters)
        # unit_celerity = list(map(lambda x: _global_parameters[variablesdefinition.celerite] * x / mean(unit_slope),
        #                          unit_slope))
        # compute inflow_unit as a separate method
        inflow_unit = self.get_inflow_unit()  #; % Débit d'entrée latéral (m3/s)
        try:
            code_transfert = _global_parameters[variablesdefinition.code_transfert]
        except KeyError:
            code_transfert = 1

        #for i in range(int(_local_parameters[variablesdefinition.nb_unit])):  # % garde fou : verifier Nb_unit > 1
        celerite_bief = _global_parameters[variablesdefinition.celerite]
        hayami_core = routing.hayamimodel(_global_parameters, unit_length, celerite_bief)
        # % Calcul du débit entrant dans un tronçon : entrée latérale (ruissellement interrill) + sortie du tronçon amont
        #if i == 0:
        inflow = list(map(lambda x: 0.5 * x, inflow_unit[variablesdefinition.streamflow_label_custom].values))
        # else:
        #     inflow = list(np.add(list(inflow_unit[variablesdefinition.streamflow_label_custom].values),
        #                          outflow_unit))

        if code_transfert == 1:
            print("2. MHYDAS_UH : Fonction de transfert par le modèle de l'onde diffusante")
            outflow_unit = routing.hayamitransfer(inflow, hayami_core, _global_parameters[variablesdefinition.dt])
        # % ********************* Crank Nicholson (numérique) ***********************
        elif code_transfert == 2:
            print("2. MHYDAS_UH : Fonction de transfert par le schéma de Crank-Nicholson")
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

        #print("Erosion 2", i, max(inflow))
        H_entree = []
        V_entree = []
        # % Calcul des hauteurs (m) et vitesses (m/s) dans la rigole en entrée d'unité élémentaire
        for j in range(len(inflow)):
            if inflow[j] == 0:
                H_entree.append(0)
                V_entree.append(0)
            else:
                rill_width = _local_parameters[variablesdefinition.larg_rill]
                #print("erosion", inflow[j], _local_parameters[variablesdefinition.rugo_strickler],
                    #rill_width, unit_slope[i], _local_parameters[variablesdefinition.seuil_conv])
                _h_entree = rillproperties.water_depth(
                    inflow[j], _local_parameters[variablesdefinition.rugo_strickler],
                    rill_width, unit_slope, _local_parameters[variablesdefinition.seuil_conv]
                )

                rill_section = rill_width * _h_entree
                _v_entree = rillproperties.water_speed(inflow[j], rill_section)

                H_entree.append(_h_entree)
                V_entree.append(_v_entree)

        #print("Erosion 3", i, max(outflow_unit))
        # % Débit (m3/s) entrant dans un tronçon : entrée latérale (ruissellement interrill) + sortie du tronçon amont
        #if i == 0:

        Q_CALC_UNIT = np.c_[np.array(inflow)]
        # % Hauteur (m) d'eau en entrée du tronçon
        H_CALC_UNIT = np.c_[np.array(H_entree)]
        # % Vitesse (m/s) en entrée du tronçon
        V_CALC_UNIT = np.c_[np.array(V_entree)]
        # % Débit (m3/s) sortant du tronçon considéré (débit d'entrée, Q_entree, transféré)
        Q_CALC_SORTIE  = np.c_[np.array(outflow_unit)]

        # else:
        #     Q_CALC_UNIT = np.c_[Q_CALC_UNIT, np.array(inflow)]
        #     # % Hauteur (m) d'eau en entrée du tronçon
        #     H_CALC_UNIT = np.c_[H_CALC_UNIT, np.array(H_entree)]
        #     # % Vitesse (m/s) en entrée du tronçon
        #     V_CALC_UNIT = np.c_[V_CALC_UNIT, np.array(V_entree)]
        #     # % Débit (m3/s) sortant du tronçon considéré (débit d'entrée, Q_entree, transféré)
        #     Q_CALC_SORTIE = np.c_[Q_CALC_SORTIE, np.array(outflow_unit)]

            # % ********** Calcul des variables HYDRO (communs aux 2 méthodes) **********
            # % Ecriture matrices
        # % Débit (m3/s) sortant d'un motif
        #print("Outflow", max(outflow_unit), max(inflow_unit[variablesdefinition.streamflow_label_custom].values))
        self.Q_sortie_motif = np.add(outflow_unit,
                                     list(map(lambda flow: 0.5 * flow,
                                              inflow_unit[variablesdefinition.streamflow_label_custom].values
                                              )
                                          )
                                     )
        self.Q_sortie_parcelle = list(map(lambda value: value * _local_parameters[variablesdefinition.nb_motifs],
                                     self.Q_sortie_motif))

        #print("Erosion 2", sum(self.Q_sortie_parcelle), max(self.Q_sortie_parcelle), max(self.Q_sortie_motif))
        # % Calcul des hauteurs et vitesses dans la rigole en sortie d'un motif
        H_sortie_motif = []
        V_sortie_motif = []
        for j in range(len(inflow)):
            if self.Q_sortie_motif[j] == 0:
                H_sortie_motif.insert(j, 0)
                V_sortie_motif.insert(j, 0)
            else:
                rill_width = _local_parameters[variablesdefinition.larg_rill]
                _water_depth = rillproperties.water_depth(self.Q_sortie_motif[j],
                                                          _local_parameters[variablesdefinition.rugo_strickler],
                                                          rill_width,
                                                          unit_slope,
                                                          _local_parameters[variablesdefinition.seuil_conv]
                                                          )
                H_sortie_motif.append(_water_depth)
                _section = _water_depth * rill_width
                V_sortie_motif.append(rillproperties.water_speed(self.Q_sortie_motif[j], _section))

        # % Ajout dans les  matrices (colonne Nb_unit+1) des valeurs de Q,H,V, pour le dernier noeud
        Q_CALC_UNIT = np.c_[Q_CALC_UNIT, np.array(self.Q_sortie_motif)]
        #Q_CALC_UNIT.insert(int(_local_parameters[variablesdefinition.nb_unit]), self.Q_sortie_motif)
        H_CALC_UNIT = np.c_[H_CALC_UNIT, np.array(H_sortie_motif)] # ; % Calcul avec pente du dernier tronçon
        V_CALC_UNIT = np.c_[V_CALC_UNIT, np.array(V_sortie_motif)]

        return Q_CALC_UNIT, H_CALC_UNIT, V_CALC_UNIT, Q_CALC_SORTIE

    def set_sediment_production_parameters(self):
        # if erosion == 1 :
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
        # kteste = 1
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        # unit_length = (_local_parameters[variablesdefinition.dist_1] + _local_parameters[variablesdefinition.dist_2] +
        #                _local_parameters[variablesdefinition.dist_3]) / _local_parameters[
        #                   variablesdefinition.nb_unit]  # ; % Longueur de l'unité élémentaire
        unit_length = _local_parameters[variablesdefinition.long_uh]                            
        # unit_width = _local_parameters[variablesdefinition.larg_uh] / _local_parameters[
        #     variablesdefinition.nb_motifs]  # ; % Largeur de l'unité élémentaire (= largeur motif élémentaire)
        unit_area = _local_parameters[variablesdefinition.aire_uh] / _local_parameters[
             variablesdefinition.nb_motifs]        
        unit_slope = _local_parameters[variablesdefinition.pente_uh]#slopes.vectorize(_local_parameters)

        #precipitation_data = self.get_precipitation_data()
        ke_nu = kineticenergy.rain_on_bare_soil(self.precipitation_data, _global_parameters, kinetic_method)
        ke_veg = kineticenergy.rain_on_vegetation(_local_parameters)
        SDR_ajust_Q = sediments.adjusted_sediment_discharge_ratio(_global_parameters, _local_parameters,
                                                                  self.net_precipitation)

        self.Q_CALC_UNIT, self.H_CALC_UNIT, self.V_CALC_UNIT, Q_CALC_SORTIE = self.get_inflow_from_net_precipitation()

        type_splash = 2
        #precipitation_data = self.get_precipitation_data()
        inflow_unit = self.get_inflow_unit()
        if type_splash == 1:
            self.splash_method = 'LISEM'
            print('3. MHYDAS_UH_EROSION : LISEM + SDR  ')
            #self.SPLASH_CALC_UNIT_LISEM, self.CONC_SPLASH, 
            #self.Splash_Unit_LISEM, self.Splash_direct, 
            #self.Splash_indirect = f_MHYDAS_UH_Splash_LISEM_Erosion(precipitation_data, Q_entree_unit, Q_CALC_UNIT, Larg_unit, Long_unit, SDR_ajust_Q, PARAM,KE_nu, KE_couv, kteste)
        if type_splash == 2:
            self.splash_method = 'MWD + SDR'
            print('3. MHYDAS_UH_EROSION : Calcul du SPLASH par MWD avec ajout SDR')
            self.SPLASH_CALC_UNIT_LISEM, self.CONC_SPLASH, self.Splash_Unit_LISEM, self.Splash_direct, self.Splash_indirect = splashmodels.mean_weight_diameter(self.precipitation_data, inflow_unit, self.Q_CALC_UNIT,
                                                                     unit_area, SDR_ajust_Q,
                                                                     _global_parameters, _local_parameters,
                                                                     ke_nu, ke_veg, self.kteste)
        elif type_splash == 3:
             self.splash_method = 'EUROSEM'
             print('3. MHYDAS_UH_EROSION : Calcul du SPLASH par EUROSEM ')
        #     SPLASH_CALC_UNIT_LISEM, CONC_SPLASH,  Splash_Unit_LISEM, Splash_direct, Splash_indirect = f_MHYDAS_UH_Splash_EUROSEM_Erosion(Pluie, Q_entree_unit, Q_CALC_UNIT, Larg_unit, Long_unit, SDR_ajust_Q, PARAM, pn,KE_nu, KE_couv)
        #        #%splash_mesure (ki, length(Pluie.int)) = (Splash_Unit_LISEM);

        # %disp ('3bis. MHYDAS_UH_EROSION : Calcul du SPLASH par WEPP');
        # %[SPLASH_CALC_UNIT_WEPP] = f_MHYDAS_UH_Splash_WEPP_Erosion(Pluie, PARAM.Argile, PARAM.Sable, PARAM.Sablefin, Larg_unit, Long_unit, pn, pente_unit, SDR_pente,PARAM.dt, PARAM);

        # %******************    CALCUL CHARGE SOLIDE DANS TRONCONS   ***************

        print('4. MHYDAS_UH_EROSION : Calcul des concentrations en MES dans les tronçons')  # ;

        self.CALC_CONC_TR_LISEM, self.CALC_PROD_INTERNE, self.MASSE_SED, self.CALC_TC_LISEM, self.CALC_VOL_TR_LISEM = \
            sediments.sediments_concentration_per_unit(self.SPLASH_CALC_UNIT_LISEM, self.H_CALC_UNIT, self.V_CALC_UNIT,
                                                       self.Q_CALC_UNIT, Q_CALC_SORTIE, unit_area, unit_length,
                                                       unit_slope, self.net_precipitation,
                                                       _global_parameters, _local_parameters
                                                       )

    def get_hydrologic_balance(self):
        # %*************************   BILAN HYDROLOGIQUE   *************************
        # % Calcul des volumes : Pluie totale, Lame infiltrée, lame ruisselée, Volume mesuré, volume calculé
        _global_parameters = self.get_global_parameters()
        #_local_parameters = self.get_local_parameters()
        #precipitation_data = self.get_precipitation_data()
        if not len(self.net_precipitation):
            self.get_infiltration_data()
        self.L_Pluie = sum(self.precipitation_data[variablesdefinition.precipitation_label_custom].values) * 1000
        self.L_Ruiss = sum(self.net_precipitation[variablesdefinition.precipitation_label_custom]) * 1000
        self.L_Inf = self.L_Pluie - self.L_Ruiss  # ; % Lames précipitées, ruisselées et infiltrées en mm
        streamflow = self.streamflow_data

        self.Vol_mes = np.trapz(y=streamflow[variablesdefinition.streamflow_label_custom].values, dx=60)  # ; % Calcul du volume écoulé mesuré en m3

        self.Vol_cal = sum(self.Q_sortie_parcelle) * _global_parameters[variablesdefinition.dt]
        _, self.Qmax_mes = max(enumerate(streamflow[variablesdefinition.streamflow_label_custom].values),
                             key=operator.itemgetter(1))
        _, self.Qmax_cal = max(enumerate(self.Q_sortie_parcelle), key=operator.itemgetter(1))

        aux = set(self.precipitation_data[variablesdefinition.datetime].values).intersection(
            streamflow[variablesdefinition.datetime].values)

        ii = self.precipitation_data.loc[self.precipitation_data[variablesdefinition.datetime].isin(aux)].index.values
        jj = streamflow.loc[streamflow[variablesdefinition.datetime].isin(aux)].index.values
        self.coeff_Nash, self.RMSE, self.CRM = errors.nash_criteria([streamflow[variablesdefinition.streamflow_label_custom].values[j]
                                                      for j in jj],
                                                      [self.Q_sortie_parcelle[i] for i in ii])

    def set_erosion_balance_parameters(self):
        # %*******************   BILAN EROSION (mesurée & simulée) *******************
        # % Bilan au tronçon & à la parcelle
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        #streamflow = self.streamflow_data
        self.set_sediment_production_parameters()

        self.CALC_Prod_interne_Tr, self.CALC_Sortie_MES_Parcelle, self.CALC_Splash_Direct_Tot_Parcelle, \
        self.CALC_Splash_Indirect_Tot_Parcelle, self.CALC_Splash_Effectif_Parcelle = \
            sediments.sediments_balance(self.Q_CALC_UNIT,
                                        self.CALC_CONC_TR_LISEM,
                                        self.CALC_PROD_INTERNE,
                                        self.Splash_direct,
                                        self.Splash_indirect,
                                        self.SPLASH_CALC_UNIT_LISEM,
                                        _local_parameters,
                                        _global_parameters)
        # % Calcul des concentrations max (mesurée et simulée)
        self.Cmax_mes = max(self.sediment_concentration_data[variablesdefinition.concentration_label].values)
        self.Cmax_cal = max(self.CALC_CONC_TR_LISEM[:, int(_local_parameters[variablesdefinition.nb_unit])])
        # %calcul de la masse en sédiment mesurée exportée
        self.sed_mes = sediments.measured_sediment_mass(self.streamflow_data, self.sediment_concentration_data)

    def create_sedimentograph(self):
        _global_parameters = self.get_global_parameters()
        _local_parameters = self.get_local_parameters()
        #precipitation_data = self.get_precipitation_data()
        _, _, _, infiltration_data = self.get_infiltration_data()
        #streamflow_data = self.streamflow_data

        self.set_erosion_balance_parameters()
        self.get_hydrologic_balance()

        data = graphics.sedimentograph(self.precipitation_data, infiltration_data, 
                                self.streamflow_data, self.Q_sortie_parcelle, 
                                self.sediment_concentration_data, self.sed_mes, 
                                self.CALC_CONC_TR_LISEM, 
                                self.CALC_Sortie_MES_Parcelle,
                                self.CALC_Prod_interne_Tr, self.Cmax_mes, 
                                self.Cmax_cal, self.L_Pluie, self.L_Inf, 
                                self.L_Ruiss, self.Vol_mes, self.Vol_cal, 
                                self.Qmax_mes, self.Qmax_cal, self.coeff_Nash,
                                _global_parameters, _local_parameters
                                )
        graphics.erosion_balance_per_block(self.splash_method, 
                                           self.CALC_Prod_interne_Tr, 
                                           _local_parameters,
                                           _global_parameters, self.sed_mes, 
                                           self.CALC_Sortie_MES_Parcelle,
                                           self.CALC_Splash_Effectif_Parcelle, 
                                           self.CALC_Splash_Direct_Tot_Parcelle,
                                           self.CALC_Splash_Indirect_Tot_Parcelle,
                                           self.L_Pluie, self.L_Inf, 
                                           self.Vol_mes, self.Vol_cal, 
                                           self.Qmax_mes, self.Qmax_cal, 
                                           self.L_Ruiss, self.coeff_Nash
                                           )
        return data

# if __name__ == '__main__':
#     new_model = Model(data_directory=default_data_file_dir)
#     new_model.set_parameters()
#     #new_model.set_sediment_production_parameters()
#     #new_model.get_hydrologic_balance()
#     #new_model.get_inflow_from_net_precipitation()
#     #new_model.get_hydrologic_balance()
#     new_model.create_sedimentograph()
    # print(new_model.get_precipitation_data())
    # print(new_model.get_streamflow_data())
    # print(new_model.get_sediment_concentration_data())
    # new_model.get_infiltration_data()
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
