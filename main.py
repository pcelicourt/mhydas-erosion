
import glob
import os
from mhydas.mhydas.erosion import erosion


data_file_dir = "./mhydas/mhydas/data/"

new_erosion_model = erosion.Model(data_file_dir)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    config_file = glob.glob(data_file_dir + "config.ini")
    if os.path.exists(config_file[-1]):
        print(glob.glob(data_file_dir + "config.ini"))
    new_erosion_model.set_parameters()
    new_erosion_model.create_sedimentograph()
