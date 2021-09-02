import pandas as pd
import numpy as np
import math
import operator

from mhydas.mhydas.utilities import variablesdefinition


def morelseytouxmethod(precipitation, storage_and_suction_factor, parameters_as_dict):
    if isinstance(precipitation, pd.DataFrame):
        distance_resolution = 0.00000005
        precipitation_data_length = precipitation.shape[0]
        infiltration_capacity = np.zeros(precipitation_data_length)
        _net_precipitation = np.zeros(precipitation_data_length)
        r = list(map(lambda _value: _value / parameters_as_dict[variablesdefinition.dt],
                precipitation[variablesdefinition.precipitation_label_custom].values))
        re = list(map(lambda _value: _value / parameters_as_dict[variablesdefinition.ks], r))

        cumulative_precipitation = 0
        cond_sat = 0
        _precipitation_values = precipitation[variablesdefinition.precipitation_label_custom].values

        index = 0
        while cond_sat < 1 and index <= (precipitation_data_length - 1):
            index += 1
            cumulative_precipitation += _precipitation_values[index-1]
            if re[index] > 1: #tp1 = (ii - 1) * PARAM.dt + (1 / r(ii)) * ((Sf / (re(ii) - 1) - som))
                tp1 = index * parameters_as_dict[variablesdefinition.dt] + (1 / r[index]) * \
                      (storage_and_suction_factor / (re[index] - 1) - cumulative_precipitation)
                if tp1 <= (index+1) * parameters_as_dict[variablesdefinition.dt]: #(1 or 2)
                    if tp1 > (index) * parameters_as_dict[variablesdefinition.dt]:# (index+1)
                        cas = 1
                        ip = index
                        tp = tp1
                        rp = r[index]
                        wp = cumulative_precipitation + \
                             (tp1 - (index) * parameters_as_dict[variablesdefinition.dt]) * r[index]
                        cond_sat = 2
                    else:
                        cas = 2
                        tp = (index) * parameters_as_dict[variablesdefinition.dt]
                        ip = index-1
                        rp = (1 + storage_and_suction_factor / cumulative_precipitation) * \
                             parameters_as_dict[variablesdefinition.ks]
                        wp = cumulative_precipitation
                        cond_sat = 3
        # while cond_sat < 1 and index <= (precipitation_data_length - 1):
        #     index += 1
        #     cumulative_precipitation += _precipitation_values[index-2]
        #     if re[index-1] > 1: #tp1 = (ii - 1) * PARAM.dt + (1 / r(ii)) * ((Sf / (re(ii) - 1) - som))
        #         tp1 = (index-1) * parameters_as_dict[variablesdefinition.dt] + (1 / r[index-1]) * \
        #               (storage_and_suction_factor / (re[index-1] - 1) - cumulative_precipitation)
        #         if tp1 <= index * parameters_as_dict[variablesdefinition.dt]: #(1 or 2)
        #             if tp1 > (index-1) * parameters_as_dict[variablesdefinition.dt]:# (index+1)
        #                 cas = 1
        #                 ip = index
        #                 tp = tp1
        #                 rp = r[index-1]
        #                 wp = cumulative_precipitation + \
        #                      (tp1 - (index-1) * parameters_as_dict[variablesdefinition.dt]) * r[index-1]
        #                 cond_sat = 2
        #             else:
        #                 cas = 2
        #                 tp = (index-1) * parameters_as_dict[variablesdefinition.dt]
        #                 ip = index-1
        #                 rp = (1 + storage_and_suction_factor / cumulative_precipitation) * \
        #                      parameters_as_dict[variablesdefinition.ks]
        #                 wp = cumulative_precipitation
        #                 cond_sat = 3
        time_steps = list(map(lambda i: (i+1) * parameters_as_dict[variablesdefinition.dt],
                                  range(precipitation_data_length)))

        #dw = np.zeros(precipitation_data_length)

        if cond_sat == 0:
            #for i in range(precipitation_data_length):
            # time_steps = list(map(lambda i: i * parameters_as_dict[variablesdefinition.dt],
            #                       range(precipitation_data_length)))
            _net_precipitation = np.zeros(precipitation_data_length)
            infiltration_capacity = np.zeros(precipitation_data_length)
        else:
            rpe = rp / parameters_as_dict[variablesdefinition.ks]
            deltawi = 0
            dw1 = 0

            #dw = []
            innertime_steps = [0]*precipitation_data_length
            dw = [0]*precipitation_data_length
            for i in range(ip+1, precipitation_data_length):
                innertime_steps.insert(i, (i+1) * parameters_as_dict[variablesdefinition.dt] - tp)
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
                    #print("infiltration", i, innertime_steps[i], parameters_as_dict["dt2"], parameters_as_dict["dt1"])
                    if innertime_steps[i]<=parameters_as_dict["dt2"] and innertime_steps[i]>parameters_as_dict["dt1"]:
                        deltawi = dw2 - 0.5 * distance_resolution
                        dw1 = dw2 - 2 * distance_resolution
                        critere = 1
                    dw2 += distance_resolution
                    parameters_as_dict["dt1"] = parameters_as_dict["dt2"]

                dw.insert(i, deltawi)
                infiltration_capacity[i] = (dw[i] - dw[i-1]) / (innertime_steps[i] - innertime_steps[i-1])

            _net_precipitation = [0]*(ip+1)
            for i in range(ip+1, precipitation_data_length):
                if r[i] > infiltration_capacity[i]:
                    _net_precipitation.append((r[i] - infiltration_capacity[i]) *
                                              parameters_as_dict[variablesdefinition.dt])
                if r[i] <= infiltration_capacity[i]:
                    _net_precipitation.append(0)
        infiltration_capacity = pd.DataFrame({variablesdefinition.datetime: precipitation[variablesdefinition.datetime].values,
                                     variablesdefinition.infiltration_rate_label_custom: infiltration_capacity})
        net_precipitation = pd.DataFrame({variablesdefinition.datetime:
                                              precipitation[variablesdefinition.datetime].values,
                                     variablesdefinition.precipitation_label_custom: _net_precipitation})
        infiltration = pd.DataFrame({variablesdefinition.datetime:
                                              precipitation[variablesdefinition.datetime].values,
                                     variablesdefinition.infiltration_rate_label_custom: list(map(operator.sub,
                                                                                                  _precipitation_values,
                                                                                                  _net_precipitation))})

        return time_steps, infiltration_capacity, net_precipitation, infiltration

def philipmethod():
    pass
