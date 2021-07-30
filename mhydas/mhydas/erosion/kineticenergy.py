import pandas as pd
import numpy as np
import math

from mhydas.mhydas.utilities import variablesdefinition


def rain_on_bare_soil(precipitation, parameters_as_dict, kinetic_method):
    #function [KE_nu,KE_couv] = f_MHYDAS_UH_Kinetic_Energy(Pluie,PARAM,kinetic_method)

    # % Autor: S.Gumiere
    # % data : 12/02/2007
    # %
    # % Descrição: Calculo da Energia cinetica das gotas de chuva por 4 metodos da literatura
    # % utilizado como base a tese de Tomas(1992)e os trabalhos de:
    # % Hudson(1986)-Soil concervation, London
    # % Kinnell(1981)- Rainfall intensity-Kinetic energy relationships for soil
    # % loss pretiction - SSSJA - v.45
    # % Rosewell(1986) - Rainfall Kinetic energy in Eastern Australia - American
    # % Meteorol. Society.
    # %Salles&Poesen (2000)-Rain properties controlling sil splash detachment HP
    #
    # %Entrada : Pluie,PARAM,kinetic_method
    # %Saida   : KE --> Kinetic Energy (J*m^2*mm^-1)
    ke_bs = []
    kinetic_energy_methods = {
        1: lambda x: 9.81 + 11.25*math.log10(x),
        2: lambda x: 29.8 * (1-(4.29/x)),
        3: lambda x: 28.12 if x > 76.0 else 11.9 + 8.73*math.log10(x),
        4: lambda x: 29.22*(1-0.894*np.exp(-0.0477*x)),
        5: lambda x: 29.31*(1-0.281*np.exp(-0.018*x)),
        6: lambda x: 8.95 + 8.44*math.log10(x),
        7: lambda x: 8.95 + 0.5546*x- 0.5009e-2*x**2 + 0.126e-4*x**3,
        8: lambda x: 27.30 + 21.68*np.exp(-0.048*x)-41.26*np.exp(-0.072*x),
        9: lambda x: 26.35*(1-0.669*np.exp(-0.0349*x)),
        10: lambda x: 29*(1-0.596*np.exp(-0.0404*x)),
        11: lambda x: 24.48 * (1-1.253/x),
        12: lambda x: 24.8 * (1-1.292/x),
        13: lambda x: 29.86 * (1-4.288/x),
        14: lambda x: 30.13 * (1-5.482/x)
    }
    if isinstance(precipitation, pd.DataFrame):
        for i in range(parameters_as_dict[variablesdefinition.nb_dt]):
            for j in range(len(precipitation)):
                if precipitation[variablesdefinition.precipitation_label].values[j] < 0.0000001:
                    ke_bs.append(0)
                else:
                    ke_bs.append(kinetic_energy_methods[kinetic_method](
                       precipitation[variablesdefinition.precipitation_label].values[j]
                       )
                    )
        #ke_couv = 15.8*(parameters_as_dict[variablesdefinition.haut_canopee])**0.5-5.87
        return ke_bs
    else:
        raise ValueError("The precipitation data is not a data frame.")


def rain_on_vegetation(parameters_as_dict):
    ke_veg = 15.8 * math.sqrt(parameters_as_dict[variablesdefinition.haut_canopee]) - 5.87
    return ke_veg
