import pandas as pd
from datetime import datetime

from mhydas.mhydas.utilities import variablesdefinition

number_of_seconds_per_minute = 60


#This function only works for a csv file formatted as
#%Annee Mois Jour Heure Minute Sec ROG3
#% - - - - - - (mm/h) (mm/h) (mm/h)
#1997	10	7	13	45	0	 0.50
def disaggregate_date_time_from_minute_to_seconds(data, dt):
    columns = [variablesdefinition.datetime, variablesdefinition.precipitation_label,
               variablesdefinition.precipitation_label_custom
            ]
    precipitation_data = pd.DataFrame(columns=columns)

    if number_of_seconds_per_minute % dt == 0:
        number_of_aggregates_per_minute = int(number_of_seconds_per_minute/dt)
        #data = pd.read_csv(file_path, sep=delimiter, skiprows=2, header=None).convert_dtypes()
        for index, row in data.iterrows():
            for number_of_row in range(number_of_aggregates_per_minute):
                _time = datetime(*(row.values[:5]), int(number_of_row * dt))
                values = [_time, row.values[6], row.values[6] * 0.001 * dt/3600
                          ]
                
                #precipitation_data = precipitation_data.append(dict(zip(columns, values)), ignore_index=True)
                precipitation_data = pd.concat([precipitation_data, 
                                                pd.DataFrame.from_records([dict(zip(columns, values))])], 
                                               ignore_index=True, sort=False)
                
        return precipitation_data
    else:
        raise ValueError("The time step parameter value is not a divisor of 60.")
