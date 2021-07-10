import glob
import os
import errno
from datetime import datetime

import pandas as pd

from mhydas.mhydas.utilities import filesmanager, variablesdefinition
from mhydas.mhydas.utilities.datetimetransformation import disaggregate_date_time_from_minute_to_seconds

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
        columns = ['Datetime', 'flow_L_s', 'flow_m3_s']
        streamflow_data = pd.DataFrame(columns=columns)
        data = pd.read_csv(self.main_config_file_content.get(variablesdefinition.streamflow),
                                        sep="\t", skiprows=2
                           )
        for index, row in data.iterrows():
            values = [datetime(*(row.values[:6])), row.values[6], row.values[6] * 0.001]
            streamflow_data = streamflow_data.append(dict(zip(columns, values)), ignore_index=True)
        return streamflow_data

if __name__ == '__main__':
    new_model = Model(data_directory=default_data_file_dir)
    new_model.set_parameters()
    #print(new_model.get_precipitation_data())
    print(new_model.get_streamflow_data())
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