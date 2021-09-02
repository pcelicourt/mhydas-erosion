import math
import numpy as np

def nash_criteria(Q_mes, Q_cal):
    # function [Critere_Nash, RMSE, CRM] = f_MHYDAS_UH_Nash(Q_mes,Q_cal)
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Calcul du critère de Nash
    # % Nash JE, Sutcliffe JV. 1970. River flow forecasting through conceptual model. Part I: a discussion of principles. Journal of Hydrology 10: 282-290.
    #
    # %*********************    SORTIES DE LA FONCTION    ***********************
    # % Critère_Nash : quantifie l'erreur pour une simulation donnée --> différence entre les courbes simulée & observée
    #
    # %********************   ARGUMENTS DE LA FONCTION    ***********************
    # % Q_mes : débit mesuré
    # % Q_cal : débit calculé

    Nb_Q = len(Q_mes)
    mean_measured_flow = sum(Q_mes) / Nb_Q

    num = sum(list(map(lambda dif: dif**2, np.subtract(Q_mes, Q_cal))))
    den = sum(list(map(lambda measured_flow: (measured_flow-mean_measured_flow)**2, Q_mes)))

    nash_criteria = 1-(num / den)

    #%**************************************************************************
    #% Calculo do RMSE --> Silvio

    RMSE = math.sqrt(num / Nb_Q)*100

    #%**************************************************************************
    #% Calculo coeficiente de Massa Residual --> Silvio

    crm = (sum(Q_mes) - sum(Q_cal))/sum(Q_mes)

    #%**************************************************************************
    #% Calculo do erro maximo --> Silvio

    #ME = max(abs(Q_mes - Q_cal))
    return nash_criteria, RMSE, crm
    #%**************************************************************************