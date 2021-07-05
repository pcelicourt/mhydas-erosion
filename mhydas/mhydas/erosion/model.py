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
        self.default_global_param_config_file = "global_prameters.ini"
        self.default_specific_param_config_file = "specific_prameters.ini"

    def check_file_existence(self, file_path):
        if os.path.exists(file_path):
            return file_path
        else:
            return None

    def set_config_file(self, file_path, default_path=None):
        if file_path:
            return self.check_file_existence(file_path)
        else:
            default_config_file_path = glob.glob(default_data_file_dir + self.default_config_file)[-1]
            return self.check_file_existence(default_config_file_path)

    def set_main_config_file(self, file_path=None):
        if file_path:
            self.main_config_file_path = self.set_config_file(file_path)
        else:
            default_config_file_path = glob.glob(default_data_file_dir + self.default_main_config_file)[-1]
            self.main_config_file_path = self.set_config_file(default_config_file_path)

    def set_global_parameters_config_file(self, file_path=None):
        if file_path:
            self.global_param_config_file_path = self.set_config_file(file_path)
        else:
            default_config_file_path = glob.glob(default_data_file_dir + self.default_global_param_config_file)[-1]
            self.global_param_config_file_path = self.set_config_file(default_config_file_path)

    def set_specific_parameters_config_file(self, param_file_path=None):
        if config_file_path:
            self.specific_param_config_file_path = self.set_config_file(config_file_path)
        else:
            default_config_file_path = glob.glob(default_data_file_dir + self.default_specific_param_config_file)[-1]
            self.specific_param_config_file_path = self.set_config_file(default_config_file_path)

    def get_config_file(self):
        return self.config_file_path

    def get_global_parameters_file(self):
        return self.global_param_file_path

    def get_specific_parameters_config_file(self):
        return self.specific_param_file_path

if __name__ == '__main__':
    config_file = glob.glob(default_data_file_dir + "/*.ini")
    new_model = Model(data_directory=default_data_file_dir)
    new_model.set_main_config_file()
    if os.path.exists(config_file[0]):
        print(filesmanager.read_main_config_file(new_model.main_config_file_path))