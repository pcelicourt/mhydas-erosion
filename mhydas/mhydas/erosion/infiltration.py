import pandas as pd
import numpy as np
import math

from mhydas.mhydas.utilities import variablesdefinition

def morelseytouxmethod(precipitation, storage_and_suction_factor, parameters_as_dict):
    if isinstance(precipitation, pd.DataFrame):
        distance_resolution = 0.00000005
        precipitation_data_length = precipitation.shape[0]
        infiltration_capacity = np.zeros(precipitation_data_length)
        net_precipitation = np.zeros(precipitation_data_length)
        r = precipitation[variablesdefinition.precipitation_label_custom] / parameters_as_dict[variablesdefinition.dt]
        re = r / parameters_as_dict[variablesdefinition.ks]
        cumulative_precipitation = 0
        cond_sat = 0
        for index, value in precipitation.iterrows():
            if cond_sat < 1 and index <= (precipitation_data_length - 1):
                cumulative_precipitation += precipitation[variablesdefinition.precipitation_label_custom].values[index]
                if re[index + 1] > 1:
                    tp1 = index * parameters_as_dict[variablesdefinition.dt] + (1 / r[index + 1]) * \
                          ((storage_and_suction_factor / (re[index + 1] - 1)) - cumulative_precipitation)
                    if tp1 <= (index + 1) * parameters_as_dict[variablesdefinition.dt]:
                        if tp1 > index * parameters_as_dict[variablesdefinition.dt]:
                            cas = 1
                            ip = index + 1
                            tp = tp1
                            rp = r[index + 1]
                            wp = cumulative_precipitation + \
                                 (tp1 - index * parameters_as_dict[variablesdefinition.dt]) * r[index + 1]
                            cond_sat = 2
                        else:
                            cas = 2
                            tp = index * parameters_as_dict[variablesdefinition.dt]
                            ip = index
                            rp = (1 + (storage_and_suction_factor / cumulative_precipitation)) * \
                                 parameters_as_dict[variablesdefinition.ks]
                            wp = cumulative_precipitation
                            cond_sat = 3
        time_steps = []
        for i in range(precipitation_data_length):
            time_steps.append(parameters_as_dict[variablesdefinition.dt] * i)
        dw = np.zeros(precipitation_data_length)
        print("cond_sat", cond_sat)
        if cond_sat == 0:
            # for i in range(precipitation_data_length):
            #     time_steps.append(i * parameters_as_dict[variablesdefinition.dt])
            net_precipitation = np.zeros(precipitation_data_length)
            infiltration_capacity = np.zeros(precipitation_data_length)
        else:
            rpe = rp / parameters_as_dict[variablesdefinition.ks]
            # for i in range(precipitation_data_length):
            #     time_steps.append(parameters_as_dict[variablesdefinition.dt] * i)
            deltawi = 0
            dw1 = 0
            innertime_steps = np.zeros(precipitation_data_length)
            for i in range(ip +1, precipitation_data_length):
                innertime_steps[i] = i * parameters_as_dict[variablesdefinition.dt] - tp
                critere = 0
                parameters_as_dict["dt1"] = dw1 * parameters_as_dict[variablesdefinition.beta] / \
                                            parameters_as_dict[variablesdefinition.ks] - \
                                            (wp * (parameters_as_dict[variablesdefinition.beta] * rpe - 1) /
                                             parameters_as_dict[variablesdefinition.ks]) * \
                                            math.log((rpe * wp + dw1) / (rpe * wp))
                dw2 = dw1 + distance_resolution
                while critere == 0:
                    parameters_as_dict["dt2"] = (dw2 * parameters_as_dict[variablesdefinition.beta] /
                                                 parameters_as_dict[variablesdefinition.ks]) - \
                                                (wp * (parameters_as_dict[variablesdefinition.beta] * rpe - 1) /
                                                 parameters_as_dict[variablesdefinition.ks]) * \
                                                math.log((rpe * wp + dw2) / (rpe * wp))
                    if innertime_steps[i]<=parameters_as_dict["dt2"] and innertime_steps[i]>parameters_as_dict["dt1"]:
                        deltawi = dw2 - 0.5 * distance_resolution
                        dw1 = dw2 - 2 * distance_resolution
                        critere = 1
                    dw2 += distance_resolution
                    parameters_as_dict["dt1"] = parameters_as_dict["dt2"]

                dw[i] = deltawi
                infiltration_capacity[i] = (dw[i] - dw[i-1]) / (innertime_steps[i] - innertime_steps[i-1])

            for i in range(ip):
                net_precipitation[i] = 0

            for i in range(ip+1, precipitation_data_length):
                if r[i] > infiltration_capacity[i]:
                    net_precipitation[i] = (r[i] - infiltration_capacity[i]) * \
                                           parameters_as_dict[variablesdefinition.dt]
                if r[i] <= infiltration_capacity[i]:
                    net_precipitation[i] = 0
        #precipitation["infiltration_capacity"] = infiltration_capacity
        #precipitation["time_step"] = time_steps
        #precipitation["net_precipitation"] = net_precipitation
        #precipitation["infiltration_rate"] = \
            #precipitation[variablesdefinition.precipitation_label_custom]-precipitation["net_precipitation"]
        return time_steps, infiltration_capacity, net_precipitation#precipitation#time_steps, infiltration_capacity, net_precipitation

def philipmethod():
    pass
