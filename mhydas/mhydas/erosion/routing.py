import numpy as np
import math

from mhydas.mhydas.utilities import variablesdefinition


def hayamimodel(parameters_as_dict, long_bief, celerite_bief):
    teta = long_bief / celerite_bief
    zed = celerite_bief * long_bief / (4 * parameters_as_dict[variablesdefinition.sigma])
    if zed < 0.5:
        zed = 0.5
    if zed > 50:
        zed = 50

    noyeau = []
    var1 = (teta * zed / 3.1416) ** 0.5
    if teta > 0.1 * parameters_as_dict[variablesdefinition.dt]:
        for i in range(1, int(parameters_as_dict[variablesdefinition.duree_max_noyau]) + 1):
            t = (i - 0.5) * parameters_as_dict[variablesdefinition.dt]
            var2 = np.exp(zed * (2 - (t / teta) - (teta / t)))
            var3 = math.pow(t, 1.5)
            noyeau.append(var1 * var2 / var3)

    else:
        noyeau[0] = 0.5 / parameters_as_dict[variablesdefinition.dt]
        noyeau[1] = 0.5 / parameters_as_dict[variablesdefinition.dt]
        for i in range(2, parameters_as_dict[variablesdefinition.duree_max_noyau]):
            noyeau.append(0)
    volume = sum(noyeau) * parameters_as_dict[variablesdefinition.dt]

    noyeau = list(map(lambda x: x / volume, noyeau))
    return noyeau
    #volume_corrige = sum(noyeau) * parameters_as_dict[variablesdefinition.dt]


def hayamitransfer(inflow, hydrologic_unit, dt):
    #  % Transfert par la méthode de l'Onde diffusante résolue par la méthode analytique d'Hayami.
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Moussa R, Bocquillon C. 1996. Criteria for the choice of flood-routing methods in natural channels. Journal of Hydrology 186(1-4) : 1-30.
    # %
    # % pn      : pluie nette (mm / pas de temps)
    # % Surface : suface du bassin en m2
    # % Q_cal   : débit de sortie (l/s)
    # % dt      : pas de temps (s) --> PARAM.dt
    #
    # Nb_pn=length(Q_entree);
    # Nb_HU=length(HU);
    #
    # %_______________________________________________
    #
    # for ii = 1 : Nb_pn
    #     Q_sortie_cal(ii,1) = 0;
    #     Nb_pas_calcul = Nb_HU;
    #     if ii < Nb_HU;
    # 	    Nb_pas_calcul = ii;
    #     end
    #     for jj=1:Nb_pas_calcul
    #         Q_sortie_cal(ii,1) = Q_sortie_cal (ii,1) + (HU(jj) * Q_entree(ii-jj+1) * PARAM.dt);
    #         if Q_sortie_cal(ii,1)<0.000001
    #            Q_sortie_cal(ii,1) = 0;
    #         end
    #    end
    # end
    #print(inflow)
    net_rainfall_length = len(inflow)
    number_of_hydrologic_unit = len(hydrologic_unit)
    outflow = []
    for i in range(net_rainfall_length):
        outflow.insert(i, 0)
        number_of_calculation_steps = number_of_hydrologic_unit
        if i < number_of_hydrologic_unit:
            number_of_calculation_steps = i
        for j in range(number_of_calculation_steps):
            outflow[i] += hydrologic_unit[j] * inflow[i-j] * dt

            if outflow[i] < 0.000001:
                outflow[i] = 0
    #print(len(outflow), net_rainfall_length)
    return outflow

def crancknicholsontransfer(inflow, celerity, sigma, surface_unit_length, time_step, rainfall_series_length,
                            number_of_element_per_surface_unit, number_of_time_steps, number_of_virtual_steps):
    # Q_sortie_unit = f_MHYDAS_UH_Transfert_CNX (inflow, celerite_bief, PARAM.sigma, Long_unit, PARAM.dt, length(Pluie.time), PARAM.nb_dx, PARAM.nb_dt, PARAM.number_of_virtual_steps);
    #function Q_cal = f_MHYDAS_UH_Transfert_CNX (Q_entree, celerite, sigma, longueur, dt_mes, n, nn1, nt1, nb_pas_fictifs)
    #% number_of_element_per_surface_unit  : Nombre de pas d'espace dx
    #% number_of_time_steps  : Nombre de pas de temps pour subdiviser le pas de temps d'entrée
    #% nn   : Nombre de pas d'esapce sur le calcul on effectue les calculs  nn = number_of_element_per_surface_unit + des pas d'espace fictifs
    #% m    : Numéro du aps d'esapce pour lequel on sort les résultats m = nn -1
    #% n    : Nombre de pas de temps de mesure (à t=0, correspond l'indice 1)
    dt_cal = time_step / number_of_time_steps
    delta_x = surface_unit_length / number_of_element_per_surface_unit
    m = number_of_element_per_surface_unit + 1
    nn = number_of_element_per_surface_unit + number_of_virtual_steps
    nt = ((rainfall_series_length - 1) * number_of_time_steps) + 1

    #inflow at t = 0
    deba = []
    for i in range(1, m+1):
        alpha = 1 - (i - 1) / (m - 1)
        beta = 1 - alpha
        deba.append((alpha * inflow[0]) + (beta * inflow[0]))
    
    for i in range(m+1,1, nn):
        deba.append(inflow[0])

    #inflow at x = 0 - -------------------------------------
    deb0 = []#deb0(1) = deb1(1);
    if number_of_time_steps == 1:
        # Cas où on ne subdivise pas le temps
        deb0 += inflow[:nt]
        #for i in range(nt):
            #deb0.append()#deb0(ii) = deb1(ii)
    else:
        #% Ca où on subdivise le aps de temps de mesure dt
       for i in range(rainfall_series_length):# % Boucle sur les aps de temps de mesures
           for j in range(number_of_time_steps):#% Boucle sur les subdivisions
               indice = i * number_of_time_steps + j#% Indice de 1 à nt
               alpha  = j/number_of_time_steps
               beta = 1 - alpha# % Coefficients de l'interpolation linéaire
               deb0[indice] = alpha * inflow[i] + beta * inflow[i+1] #% Calcul de débit à x= 0 au temps correspodant à indice

    #% ---------- Calcul des debits simules - --------------------------
    h = celerity * dt_cal / delta_x
    g = (sigma * dt_cal) / (delta_x * delta_x)
    p = (-h / 4) - (g / 2)
    q = 1 + g
    r = (h / 4) - (g / 2)
    pp = (h / 4) + (g / 2)
    qq = 1 - g
    rr = (-h / 4) + (g / 2)
    debss = []
    #debss(1) = deb2(1)
    jk = 1
    debe1 = []
    debe2 = []
    ws = [0]
    for j in range(1, nt):
        debe1[0] = deb0[j]
        for i in range(1, nn):
            debe1[i] = deba[i]
        #% calcul de x  y  et   w
        debe2[0] = debe1[0]
        #% *********** Calcul des apports ***********
        #% for ii = 2 : nn ; apport(ii,jj) = 0.0 * debe1(ii); end % Ajouter la
        #% fonction d'apports, soit une fonction proportionnelle aux débits
        #% d'entrée

        for i in range(1, nn):
            ws.append(pp*debe1[i] + qq*debe1[i+1] + rr*debe1[i+2])
            #%*********** tenir compte des apports ************
            #% ws(ii) = ws(ii) + (pp-p) * apport(ii-1,jj) + (r-rr)*apport(ii,jj);

        #% ___________________ Double balayage  de Crank Nicholson _________________
        xs = []
        ys = []
        xs.insert(nn, (-q + ((q*q - 4*p*r)**0.5)) /(2*r))#; % Condition aval
        ys.insert(nn, 0) #; % Condition aval
        #k = nn - 1

        for i in np.arange(nn - 1, -1, -1):
            #ind = nn - i
            xs.insert(i-1, -p/(q+r*xs[i]))
            ys.insert(i-1, (ws[i] - r*ys[i])/(q + r*xs[i]))

        for i in range(1, nn):
            debe2.insert(i, xs[i] * debe2[i-1] + ys[i])

        for i in range(nn):
            deba.insert(i, debe2[i])

        jjj = (math.floor((j-1)/number_of_time_steps)*number_of_time_steps) - (j-1)
        if jjj == 0:
            jk = jk + 1
            debss.insert(jk, debe2[m])
    Q_cal = np.matrix.H(debss)
    return Q_cal