import glob
import os

from mhydas.mhydas.utilities import filesmanager


default_data_file_dir = "../data/"


class Model():
    def __init__(self, data_directory=None):
        if data_directory:
            self.data_directory = data_directory
        else:
            self.data_directory = default_data_file_dir
        self.default_main_config_file = "config.ini"
        self.default_global_param_config_file = "global_parameters.ini"
        self.default_specific_param_config_file = "specific_parameters.ini"

    def check_file_existence(self, file_path):
        if os.path.exists(file_path):
            return file_path
        else:
            return None

    def set_config_file(self, file_path, default_path=None):
        if file_path:
            return self.check_file_existence(file_path)
        else:
            if default_path:
                return self.check_file_existence(default_path)
            else:
                return None

    def set_main_config_file(self, file_path=None):
        if file_path:
            self.main_config_file_path = self.set_config_file(file_path)
        else:
            default_config_file_path = glob.glob(os.path.join(self.data_directory,
                                                              self.default_main_config_file
                                                              )
                                                 )[-1]
            self.main_config_file_path = self.set_config_file(default_config_file_path)

    def set_global_parameters_config_file(self, file_path=None):
        if file_path:
            self.global_param_config_file_path = self.set_config_file(file_path)
        else:
            default_config_file_path = glob.glob(os.path.join(self.data_directory,
                                                              self.default_global_param_config_file)
                                                 )[-1]
            self.global_param_config_file_path = self.set_config_file(default_config_file_path)

    def set_specific_parameters_config_file(self, file_path=None):
        if file_path:
            self.specific_param_config_file_path = self.set_config_file(file_path)
        else:
            default_config_file_path = glob.glob(os.path.join(self.data_directory,
                                                              self.default_specific_param_config_file)
                                                 )[-1]
            self.specific_param_config_file_path = self.set_config_file(default_config_file_path)

    def get_config_file(self):
        try:
            return self.config_file_path
        except:
            raise

    def get_global_parameters_file(self):
        try:
            return self.global_param_config_file_path
        except:
            raise

    def get_specific_parameters_config_file(self):
        try:
            return self.specific_param_config_file_path
        except:
            raise

if __name__ == '__main__':
    new_model = Model(data_directory=default_data_file_dir)
    new_model.set_main_config_file()
    new_model.set_global_parameters_config_file()
    new_model.set_specific_parameters_config_file()
    print(filesmanager.read_main_config_file(new_model.main_config_file_path))
    print(filesmanager.read_model_parameters_config_file(new_model.global_param_config_file_path))
    print(filesmanager.read_model_parameters_config_file(new_model.specific_param_config_file_path))