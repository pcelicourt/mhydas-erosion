# function[t, Ic, PNette] = f_MHYDAS_UH_Morel_Seytoux(Pluie, Sf, PARAM)
#
# % Modèle
# de
# calcul
# de
# l
# 'infiltration par la méthode de Morel-Seytoux qui est une variante du modèle de Green et AMpt (1991)
# % Auteurs: Gumiere., S.J., D.Raclot & G.Davy
# % Version: 2008
# % Date: Première
# version
# en
# Fortran(1989;
# Thèse
# Moussa, 1991), Version
# Matlab(juin
# 2000).
# % Vérification
# et
# validation
# Thèse
# N.Chahinian(2004)
# % Chahinian
# N, Moussa
# R, Andrieux
# P, Voltz
# M.
# 2005.
# Comparison
# of
# infiltration
# models
# to
# simulate
# flood
# events
# at
# the
# field
# scale.Journal
# of
# Hydrology, 50
# p., in press.
# %
# % Les
# équations
# de
# base
# proviennent
# des
# articles:
# % Morel - Seytoux
# HJ.
# 1978.
# Derivation
# of
# equations
# for variable rainfall infiltration.Water Resources Research 14(4): 561 - 568.
# % Morel - Seytoux
# HJ, Khanji
# J.
# 1974.
# Derivation
# of
# an
# equation
# of
# infiltration.Water
# Resources
# Research
# 10(4): 795 - 800.
# % Morel - Seytoux
# 1980.
# Résolution
# de
# l
# 'équation de l'
# infiltration.Grenoble.
# %
# % La
# méthode
# de
# résolution
# de
# Morel - Seytoux
# a
# besoin
# d
# 'une résolution d'
# une
# équation
# par
# itération.Pour
# cela, dans
# les
# versions
# précédentes
# de
# % Matlab( < 5) on
# pourrait
# utiliser
# la
# fonciton
# fsolve.Mais
# dans
# les
# versions
# futures
# on
# n
# 'utilise plus la fonction fsolve de Matlab qui a
# % changé
# d
# 'arguments. On remplace par une résolution par un sous-programme type Newton-Raphson
#
# % PARAM.dt: pas
# de
# temps
# de
# discretisation(et
# de
# calcul) en
# secondes
# % PARAM.Ks: conductivité
# hydraulique
# à
# saturation
# en
# m / s
# % Sf = (TetaSat - TetaInitial) * H(Teta)
# en
# "m", Sf
# est
# fonction
# de
# l
# 'état initial du sol
# % PARAM.beta: Facteur
# de
# correction, constante
# entre
# 1
# et
# 1.7
# selon
# Morel - Seytoux
# et
# Khanji(1974)
# % PARAM.Ks
# et
# PARAM.beta
# sont
# les
# paramŠtres
# intrinsèques
# du
# sol, Sf
# dépend
# de
# la
# teneur
# en
# eau
# initiale
# TetaInitial
# % Pluie.int(i %): Pluie.int
# au
# pas
# i % en
# m
# % pas_resol: pas
# de
# résolution
# en
# m
# de
# l
# 'équation obtenu dans CASRES=1
#
# % __________________________________________________________________
# pas_resol = 0.00000005; % en
# m
# soit
# 0.005
# mm
# par
# pas
# de
# temps.Cette
# valeur
# devrait
# être
# mise
# comme
# paramètre
# d
# 'entrée dans le programme principal
#
# n = length(Pluie.int);
# Ic = zeros(size(Pluie.int));
# PNette = zeros(size(Pluie.int));
# r = Pluie.int. / PARAM.dt;
# re = r. / PARAM.Ks;
#
# % ** ** ** ** ** ** ** ** ** ** ** ** Calcul
# du
# temps
# de
# saturation ** ** ** ** ** ** **
# som = 0;
# cond_sat = 0;
# % Print
# "Indice temps  SomPluie.int(mm)  Temps Sat possible (h)
# ii = 1;
# while cond_sat < 1 & ii <= (n - 1)
#     ii = ii + 1;
#     som = som + Pluie.int(ii - 1);
#     if re(ii) > 1
#         tp1 = (ii - 1) * PARAM.dt + (1 / r(ii)) * ((Sf / (re(ii) - 1) - som));
#         % print
#         " i%, som, tp1", i %, som, tp1
#         if tp1 <= (ii * PARAM.dt)
#             if tp1 > ((ii - 1) * PARAM.dt)
#                 cas = 1;
#                 ip = ii;
#                 tp = tp1;
#                 rp = r(ii);
#                 wp = som + (tp1 - (ii - 1). * PARAM.dt). * r(ii);
#                 cond_sat = 2;
#             else
#                 cas = 2;
#                 tp = (ii - 1). * PARAM.dt;
#                 ip = ii - 1;
#                 rp = (1 + (Sf. / som)). * PARAM.Ks;
#                 wp = som;
#                 cond_sat = 3;
#             end
#         end
#     end
# end
#
# % ** ** ** ** ** ** ** ** ** ** Calcul
# de
# la
# capacité
# d
# 'Infiltration *****************
# if cond_sat == 0
#     for ii=1:n
#     t(ii) = ii. * PARAM.dt;
#     PNette(ii) = 0;
#     Ic(ii) = 0;
# end
#
# else
# rpe = rp. / PARAM.Ks;
# % ____________________________________________________________________
# for ii=1:n
# t(ii) = ii. * PARAM.dt;
# end
# deltawi = 0;
# dw1 = 0;
#
# for ii= ip+1:n
# temps(ii) = (ii). * PARAM.dt - tp;
# % _________________________________________________________________
# % Résolution
# de
# l
# 'équation dans le cas 1 pour MATLAB 5.3
# critere = 0;
# PARAM.dt1 = (dw1 * PARAM.beta / PARAM.Ks) - (wp * (PARAM.beta * rpe - 1) / PARAM.Ks) * log(
#     (rpe * wp + dw1) / (rpe * wp));
# dw2 = dw1 + pas_resol;
# while critere == 0
#     PARAM.dt2 = (dw2 * PARAM.beta / PARAM.Ks) - (wp * (PARAM.beta * rpe - 1) / PARAM.Ks) * log(
#         (rpe * wp + dw2) / (rpe * wp));
#     if temps(ii) <= PARAM.dt2 & temps(ii) > PARAM.dt1
#         deltawi = dw2 - 0.5 * pas_resol;
#         dw1 = dw2 - 2 * pas_resol;
#         critere = 1;
#     end
#     dw2 = dw2 + pas_resol;
#     PARAM.dt1 = PARAM.dt2;
# end
# % _________________________________________________________________
# % Fonction
# valable
# pour
# MATALAB
# 5.2
# % deltawi = fsolve('f_seytoux1', 0, [], [], PARAM.Ks, PARAM.beta, wp, rpe, temps(ii));
# % __________________________________________________________________________________
# dw(ii) = deltawi;
# Ic(ii) = (dw(ii) - dw(ii - 1)). / (temps(ii) - temps(ii - 1));
# end
#
# % ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
# % ** ** ** ** ** ** ** ** ** ** *Calcul
# de
# la
# lame
# d
# 'eau ruisselée **************
# for ii = 1:ip
# PNette(ii) = 0;
# end
# for ii = ip+1:n
# if r(ii) > Ic(ii)
#     PNette(ii) = (r(ii) - Ic(ii)) * PARAM.dt;
# end
# if r(ii) <= Ic(ii)
#     PNette(ii) = 0;
# end
# end
# end