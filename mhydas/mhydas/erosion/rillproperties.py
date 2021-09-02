import math


def water_depth(inflow, rugo_strickler, rill_width, unit_slope, seuil_conv):
    # function [hauteur_Newton,vitesse] = f_MHYDAS_UH_CalculNewtonRectangle(q,Rugo_Strickler,larg,I, seuil_conv)
    # % Calcul de la hauteur d'eau à patir du débit en applicant la méthode de
    # % Newton (convergence)
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Fichier: f_MHYDAS_UH_CalculNewton.m
    # % Les rigoles sont assimilées à des rectangles de largeur 'larg' et dont la profondeur varie avec la hauteur d'eau
    #
    # %*********************    SORTIES DE LA FONCTION    ***********************
    # % hauteur_Newton  : Hauteur d'eau (m) associée au débit q dans un tronçon calculée avec la loi de Manning-Strickler & la méthode de convergence de Newton
    # % vitesse         : Vitesse de l'eau correspondant au rapport entre le débit q et la section h*w (w : largeur de la rigole (m))
    #
    # % *******************   ARGUMENTS DE LA FONCTION    ***********************
    # % q               : débit m3/s : provient du module hydro
    # % unit_slope               : pente (m/m) : provient du vecteur pente_unit
    # % Rugo_Strickler  : rugosité strickler --> PARAM.Rugo_Strickler
    # % rill_width            : largeur de la rigole (m) --> PARAM.Larg_rill
    # % seuil_conv      : seuil de convergence (en m) --> PARAM.Seuil_Conv

    _lambda = inflow/(rugo_strickler * rill_width * math.sqrt(unit_slope)) #; % constante
    x = _lambda**(3/5) #;% valeur initiale de la hauteur d'eau
    delta = 1 + seuil_conv  #;% valeur initiale du delta (doit etre > à seuil_conv)
    #print(delta, seuil_conv, _lambda * math.pow(1/x + 2/rill_width, 2/3) - x)
    while delta > seuil_conv:
        aaa = x
        f = _lambda * (1/x + 2/rill_width)**(2/3) - x
        df = (-2/3) * _lambda * 1/x**2 * (1/x + 2/rill_width)**(-1/3)-1
        x -= f/df
        delta = abs(x - aaa)
    if x < 0.000001:
        hauteur_newton = 0
    else:
        hauteur_newton = x
        #print('rill', x)
    return hauteur_newton


def water_speed(inflow, section):
    #velocity = []
    #for i in range(len(section)):
    if section == 0:
       velocity = 0
    else:
       velocity = inflow / section
    return velocity
