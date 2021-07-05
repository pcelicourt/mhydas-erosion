
import configparser
import glob
from os import path, strerror
import errno


default_config_file_dir = "../"


def read_config_file(file_path):
    config = configparser.ConfigParser()
    file_content = {}
    if path.exists(file_path):
        model_config_file = glob.glob(file_path)
        data_directory = path.dirname(path.abspath(model_config_file[0]))
        config.read(model_config_file)
        config_file_sections = config.sections()
        for key in config_file_sections:
            for sub_key in config[key]:
                file_content[sub_key] = path.join(data_directory, config[key][sub_key])
        return file_content
    else:
        raise FileNotFoundError(
                                errno.ENOENT,
                                strerror(errno.ENOENT),
                                file_path
        )


def read_main_config_file(file_path=None):
    if not file_path:
        model_config_file = glob.glob(default_config_file_dir + "/config.ini")
    else:
        if path.exists(file_path):
            model_config_file = glob.glob(file_path)
        else:
            model_config_file = []
    return read_config_file(model_config_file[-1])


def get_global_model_parameters_config_file(file_path=None):
    if not file_path:
        model_config_file = glob.glob(default_config_file_dir + "/global_parameters.ini")
    else:
        if path.exists(file_path):
            model_config_file = glob.glob(file_path)
        else:
            model_config_file = []
    return read_config_file(model_config_file[-1])


def get_specific_model_parameters_config_file(file_path=None):
    if not file_path:
        model_config_file = glob.glob(default_config_file_dir + "/specific_parameters.ini")
    else:
        if path.exists(file_path):
            model_config_file = glob.glob(file_path)
        else:
            model_config_file = []
    return read_config_file(model_config_file[-1])