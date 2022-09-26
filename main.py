
#import glob
import os
import pandas as pd
from mhydas.mhydas.erosion import erosion, precipitation, sediments, streamflow
from mhydas.mhydas.utilities import variablesdefinition
#from mhydas.mhydas.erosion import routing, precipitation
#from mhydas.mhydas.utilities import filesmanager

data_file_dir = "./mhydas/mhydas/data/"
local_parameters_path = "local_param.csv"
global_parameters_path = "global_param.csv"

precipitation_data_path	= "Pluvio_3_1997_06_01.txt"
streamflow_data_path	= "Debit_6_1997_06_01.txt"
mes_concentration_data_path  = "Concentration_MES_NT_1997_06_01.txt"
slope_data_path = "Pente.txt"

#master_config_file = data_file_dir + "master_config.ini"
local_config_file_path = os.path.join(data_file_dir, local_parameters_path)
global_config_file_path = os.path.join(data_file_dir, global_parameters_path)
precipitation_data_path = os.path.join(data_file_dir, 
                                        precipitation_data_path)
sediment_data_path = os.path.join(data_file_dir, 
                                   mes_concentration_data_path)

streamflow_data_path = os.path.join(data_file_dir, 
                                   streamflow_data_path)
        
def execute_model():
    if os.path.exists(local_config_file_path):
        local_parameters_set = pd.read_csv(local_config_file_path).to_dict('records')
        global_parameters = pd.read_csv(global_config_file_path).to_dict('records')[0]

        precipitation_data = precipitation.disaggregate_date_time_from_minute_to_seconds(
                    precipitation_data_path,
                    "\t",
                    global_parameters[variablesdefinition.dt]
                )
        sediments_data = sediments.get_sediment_concentration_data(sediment_data_path)
        streamflow_data = streamflow.get_streamflow_data(streamflow_data_path)
        
        water_erosion_model = erosion.Model(global_parameters=global_parameters)
        water_erosion_model.set_global_parameters()
        water_erosion_model.precipitation_data = precipitation_data
        water_erosion_model.sediment_concentration_data = sediments_data
        water_erosion_model.streamflow_data = streamflow_data
        
        master_sendiment_data_frame = pd.DataFrame()
        
        for local_parameters in local_parameters_set:
            water_erosion_model.local_parameters = local_parameters
            water_erosion_model.set_local_parameters()
            sediments_data = water_erosion_model.create_sedimentograph()
            
            sediments_data = sediments_data[sediments_data["data_categories"] == "simulated erosion"]

            aggregated_sediments_data = sediments.aggregate_sediments_data(sediments_data, 
                                                                           date_column="timestamp", 
                                                                           time_step="6H"
                                                                          )
            master_sendiment_data_frame = pd.concat([master_sendiment_data_frame,
                                                     aggregated_sediments_data],
                                                    axis=1)
        master_sendiment_data_frame.columns = ["Parcelle_{0}".format(i) for i in range(1, len(local_parameters_set) + 1) ]
        sediments.save_sediments_data(master_sendiment_data_frame, "sediments.csv")
            
            #ADD RESULTS TO A CSV USING PANDAS
    else:
        water_erosion_model = erosion.Model(data_file_dir)
        
        
if __name__ == '__main__':
    execute_model()

#Press the green button in the gutter to run the script.
#new_erosion_model.set_parameters()
#new_erosion_model.create_sedimentograph()
