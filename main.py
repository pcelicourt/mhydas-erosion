import warnings
warnings.filterwarnings("ignore")

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

precipitation_data_path	= "pluvio_Philippsburg_2011.txt"
streamflow_data_path	= "debit_Philippsburg_2011.txt"
mes_concentration_data_path  = "MES_Philippsburg_2011.txt"

slope_data_path = "Pente.txt"

local_config_file_abs_path = os.path.join(data_file_dir, local_parameters_path)
global_config_file_abs_path = os.path.join(data_file_dir, global_parameters_path)


def initialize_model():
    sediment_data_abs_path = os.path.join(data_file_dir, mes_concentration_data_path)
    precipitation_data_abs_path = os.path.join(data_file_dir,precipitation_data_path)
    raw_precipitation_data = pd.read_csv(precipitation_data_abs_path, 
                                     sep="\t", skiprows=2, 
                                     header=None).convert_dtypes()
    streamflow_data_abs_path = os.path.join(data_file_dir, streamflow_data_path) 
    if os.path.exists(local_config_file_abs_path) and os.path.exists(global_config_file_abs_path):
        sediments_data = sediments.get_sediment_concentration_data(sediment_data_abs_path)
        streamflow_data = streamflow.get_streamflow_data(streamflow_data_abs_path)
        
        erosion_model = erosion.Model()
        erosion_model.raw_precipitation_data = raw_precipitation_data        
        erosion_model.sediment_concentration_data = sediments_data
        erosion_model.streamflow_data = streamflow_data
        return erosion_model
    else:
        raise 

def get_first_two_columns_of(data_frame):
    return data_frame.iloc[:, 0:2]


def get_aggregated_sediment_data(model_instance, 
                                 erosion_variable_type = variablesdefinition.simulated_erosion, 
                                 time_step="6H", show_plot=False):
    sediments_data = model_instance.create_sedimentograph(show_plot=show_plot)

    sediments_data = sediments_data[sediments_data["data_categories"] == erosion_variable_type]
    
    sediments_data = get_first_two_columns_of(sediments_data)

    aggregated_sediments_data = sediments.aggregate_sediments_data(sediments_data, 
                                                                   date_column=variablesdefinition.timestamp, 
                                                                   time_step=time_step
                                                                  )
    return aggregated_sediments_data


def execute_model(show_plot=False):
    # if os.path.exists(local_config_file_path):
    #     #TO DO: We'll need to keep track of which rows are dropped for both files
    #     local_parameters_set = pd.read_csv(local_config_file_path).to_dict('records')
    #     global_parameters = pd.read_csv(global_config_file_path).to_dict('records')#[0]
    #     #print(pd.read_csv(local_config_file_path).isnull())
    #     precipitation_data = precipitation.disaggregate_date_time_from_minute_to_seconds(
    #                 precipitation_data_path,
    #                 "\t",
    #                 global_parameters[variablesdefinition.dt]
    #             )
    #     sediments_data = sediments.get_sediment_concentration_data(sediment_data_path)
    #     streamflow_data = streamflow.get_streamflow_data(streamflow_data_path)

        erosion_model = initialize_model()

        local_parameters_set = pd.read_csv(local_config_file_abs_path).to_dict('records')
        global_parameters = pd.read_csv(global_config_file_abs_path).to_dict('records')
        
        master_sendiment_data_frame = pd.DataFrame()
        
        if len(global_parameters) == 1: 
            erosion_model.global_parameters = global_parameters[0]
            erosion_model.set_global_parameters()
            
            precipitation_data = precipitation.disaggregate_date_time_from_minute_to_seconds(
                        erosion_model.raw_precipitation_data,
                        erosion_model.global_parameters[variablesdefinition.dt]
                    )
            erosion_model.precipitation_data = precipitation_data
            
            for row_number, local_parameters in enumerate(local_parameters_set):
                erosion_model.local_parameters = local_parameters
                erosion_model.set_local_parameters()
                aggregated_sediments_data = get_aggregated_sediment_data(erosion_model)
                master_sendiment_data_frame = pd.concat([master_sendiment_data_frame,
                                                         aggregated_sediments_data],
                                                        axis=1)
        else:
            assert(len(global_parameters) == len(global_parameters))
                    
            for row_number, local_parameters in enumerate(local_parameters_set):
                erosion_model.local_parameters = local_parameters
                erosion_model.global_parameters = global_parameters[row_number]
                erosion_model.set_parameters() 
                precipitation_data = precipitation.disaggregate_date_time_from_minute_to_seconds(
                        erosion_model.raw_precipitation_data,
                        erosion_model.global_parameters[variablesdefinition.dt]
                    ) 
                erosion_model.precipitation_data = precipitation_data
                aggregated_sediments_data = get_aggregated_sediment_data(erosion_model, show_plot=show_plot)
                #print(aggregated_sediments_data)
                master_sendiment_data_frame = pd.concat([master_sendiment_data_frame,
                                                         aggregated_sediments_data],
                                                        axis=1)            
            
        #print(master_sendiment_data_frame.columns, master_sendiment_data_frame.head())
        
        master_sendiment_data_frame.columns = ["Parcelle_{0}".format(i) for i in range(1, len(local_parameters_set) + 1) ]
        sediments.save_sediments_data(master_sendiment_data_frame, "sediments.csv")
            
            #ADD RESULTS TO A CSV USING PANDAS
    #else:
        #water_erosion_model = erosion.Model(data_file_dir)
        
        
if __name__ == '__main__':
    execute_model(show_plot=False)

#Press the green button in the gutter to run the script.
#new_erosion_model.set_parameters()
#new_erosion_model.create_sedimentograph()
