# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:11:40 2022

@author: Paul Celicourt
"""
from datetime import datetime

import pandas as pd

from mhydas.mhydas.utilities import variablesdefinition

def get_streamflow_data(path):
    columns = [variablesdefinition.datetime, variablesdefinition.streamflow_label,
               variablesdefinition.streamflow_label_custom
               ]
    streamflow_data = pd.DataFrame(columns=columns)
    data = pd.read_csv(path, 
                       sep="\t", skiprows=2, header=None).convert_dtypes()
    for index, row in data.iterrows():
        _time = datetime(*(row.values[:6]))
        values = [_time, row.values[6], row.values[6] * 0.001]
        #streamflow_data = streamflow_data.append(dict(zip(columns, values)), ignore_index=True)
        streamflow_data = pd.concat([streamflow_data, 
                                                pd.DataFrame.from_records([dict(zip(columns, values))])], 
                                               ignore_index=True, sort=False)
        
    return streamflow_data