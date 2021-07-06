import configparser
import glob
import os
import errno


def read_main_config_file(file_path):
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
            print(config_file_sections)
        return file_content
    else:
        raise FileNotFoundError(
                                errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                file_path
        )


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
        raise FileNotFoundError(
                                errno.ENOENT,
                                os.strerror(errno.ENOENT),
                                file_path
        )


# def read_specggific_model_parameters_config_file(file_path=None):
#     if not file_path:
#         model_config_file = glob.glob(os.path.join(default_config_file_dir,
#                                                 "specific_parameters.ini"
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