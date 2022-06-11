
import glob
import os
from mhydas.mhydas.erosion import erosion
from mhydas.mhydas.utilities import filesmanager

data_file_dir = "./mhydas/mhydas/data/"

master_config_file = data_file_dir + "master_config.ini"


if __name__ == '__main__':
    if os.path.exists(master_config_file):
        main_config_files = filesmanager.read_master_config_file(master_config_file)
        if len(main_config_files):
            for main_config_file in main_config_files:
                configs = filesmanager.read_config_file(os.path.join(data_file_dir, main_config_file))
                new_erosion_model = erosion.Model(data_file_dir, main_config_file_name=main_config_file)
                new_erosion_model.set_parameters()
                new_erosion_model.create_sedimentograph()
    else:
        new_erosion_model = erosion.Model(data_file_dir)

#Press the green button in the gutter to run the script.
#new_erosion_model.set_parameters()
#new_erosion_model.create_sedimentograph()
