import configparser
import glob
import os
import errno


def read_master_config_file(file_path):
    config = configparser.ConfigParser()
    config_file_paths = []
    if os.path.exists(file_path):
        data_directory = os.path.dirname(os.path.abspath(file_path))
        config.read(file_path)
        config_file_sections = config.sections()
        if config_file_sections:
            for key in config_file_sections:
                for sub_key in config[key]:
                    config_files_list = config[key][sub_key]
                    sep = [x for x in [",", ";"] if x in config_files_list]
                    if len(sep):
                        sep = sep[0]
                        config_files = [_file_name.strip() for _file_name in config_files_list.split(sep)]
                        #config_file_paths = [os.path.join(data_directory, _file_name) for _file_name in config_files]
                        #print(config_file_paths)
                    else:
                        config_files = config[key][sub_key]
        else:
            pass
        return config_files
    else:
        file_not_found_error(file_path)

def read_config_file(file_path):
    config = configparser.ConfigParser()
    file_content = {}
    if os.path.exists(file_path):
        model_config_file = glob.glob(file_path)
        data_directory = os.path.dirname(os.path.abspath(model_config_file[0]))
        config.read(model_config_file)
        config_file_sections = config.sections()
        if config_file_sections:
            for key in config_file_sections:
                for sub_key in config[key]:
                    file_content[sub_key] = os.path.join(data_directory, config[key][sub_key])
        else:
            pass
        return file_content
    else:
        file_not_found_error(file_path)

def read_model_parameters_config_file(file_path=None):
    config = configparser.ConfigParser()
    file_content = {}
    if os.path.exists(file_path):
        model_config_file = glob.glob(file_path)
        data_directory = os.path.dirname(os.path.abspath(model_config_file[0]))
        config.read(model_config_file)
        config_file_sections = config.sections()
        for key in config_file_sections:
            for sub_key in config[key]:
                file_content[sub_key] = float(config[key][sub_key])
        return file_content
    else:
        file_not_found_error(file_path)

def file_not_found_error(file_name):
    raise FileNotFoundError(
        errno.ENOENT,
        os.strerror(errno.ENOENT),
        file_name
    )
# def read_specggific_model_parameters_config_file(file_path=None):
#     if not file_path:
#         model_config_file = glob.glob(os.path.join(default_config_file_dir,
#                                                 "local_parameters.ini"
#                                                 )
#                                       )
#     else:
#         if path.exists(file_path):
#             model_config_file = glob.glob(file_path)
#         else:
#             model_config_file = []
#     return read_config_file(model_config_file[-1])
#
#
#     config = configparser.ConfigParser()
#     file_content = {}
#     if os.path.exists(file_path):
#         model_config_file = glob.glob(file_path)
#         data_directory = os.path.dirname(os.path.abspath(model_config_file[0]))
#         config.read(model_config_file)
#         config_file_sections = config.sections()
#         if config_file_sections:
#             for key in config_file_sections:
#                 for sub_key in config[key]:
#                     file_content[sub_key] = os.path.join(data_directory, config[key][sub_key])
#         else:
#             print(config_file_sections)
#         return file_content
#     else:
#         raise FileNotFoundError(
#                                 errno.ENOENT,
#                                 os.strerror(errno.ENOENT),
#                                 file_path
#         )