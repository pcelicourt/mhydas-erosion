import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from mhydas.mhydas.utilities import variablesdefinition

# def hydrograph():
#     # function f_MHYDAS_UH_graphique_HYDRO(Pluie, infil, streamflow, Q_cal, PARAM)
#     # % Tracé graphique Pluie, Débit mesuré et débit calculé
#     # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
#     # % Version : 2008
#     # % Fichier: f_MHYDAS_UH_graphique_HYDRO.m
#
#     global L_Pluie   L_Inf  Vol_mes  Vol_cal   Qmax_mes  Qmax_cal  coeff_Nash
#
#     figure(1)
#     #%*************************************************************************
#     #%                                    Tracé de la pluie
#
#     subplot(2,1,1)
#     plot (Pluie.time, Pluie.int/PARAM.dt,'-b'); hold on      ; % Pluie totale
#     plot (Pluie.time,infil/PARAM.dt,'-k'); hold off ; % Infiltration
#     title (['   Pluie = ' num2str(L_Pluie) ' mm; Infiltration = ' num2str(L_Inf) ' mm']); %; Ruissellement = ' num2str(L_Ruiss) ' mm']);
#     xlabel ('Temps (HH:MM)');
#     ylabel('Pluie. et Inf. m/s');
#     legend ('Pluie','Infiltration')
#     datetick('x','HH:MM');
#     hold off
#     grid on
#
#     #%*************************************************************************
#     #%                    Tracé des hydrogrammes mesurés et calculés
#
#     subplot (2,1,2)
#     plot(Pluie.time,Q_cal,'--r');hold on; % Débit calculé en m3/s
#     plot(streamflow.time,streamflow.int,'-b') ; hold off;   % Débit mesuré en m3/s
#     title (['Vmes =' num2str(Vol_mes) ' m3; Vcal = ' num2str(Vol_cal) ' m3;    Qmax mes = ' num2str(Qmax_mes) ' m3/s;  Qmax cal = ' num2str(Qmax_cal) ' m3/s; Nash = ' num2str(coeff_Nash) ' ']);
#     xlabel ('Temps (HH:MM)');
#     ylabel ('Débit (m3/s)');
#     #%legend ('Débit calculé','Débit mesuré',1)
#     datetick('x','HH:MM');
#     #grid on
#     #hold off


def sedimentograph(Pluie, infil, streamflow, mes):
    # function f_MHYDAS_UH_graphique_MES(Pluie, infil, streamflow, MES, Q_cal, SED_mes, CALC_Sortie_MES_Parcelle, CALC_CONC_TR_LISEM, CALC_Splash_Direct_Tot_Parcelle, CALC_Splash_Indirect_Tot_Parcelle, CALC_Splash_Effectif_Parcelle, CALC_Prod_interne_Tr, PARAM,metod)
    # % Tracé graphique : - pluie, hydrogrammes & turbidigrammes mesuré et simulé (FIGURE 2)
    # %                   - tableau bilan des exportations sur la parcelle (FIGURE 3)
    # % Auteurs: Gumiere.,S.J., D. Raclot & G. Davy
    # % Version : 2008
    # % Fichier: f_MHYDAS_UH_graphique_MES.m

    #global L_Pluie  L_Inf  Vol_mes  Vol_cal  Qmax_mes  Qmax_cal Cmax_mes  Cmax_cal  L_Ruiss  coeff_Nash
    data = pd.DataFrame({"timestamp": Pluie[variablesdefinition.timestamp].values,
                       "values": Pluie[variablesdefinition.precipitation_label].values,
                       "data_group": ["precipitation"]*len(Pluie[variablesdefinition.precipitation_label]),
                       "data_categories": ["precipitation"]*len(Pluie[variablesdefinition.precipitation_label])
                       })
    data = data.append(pd.DataFrame({"timestamp": infil[variablesdefinition.timestamp].values,
                       "values": infil[variablesdefinition.infiltration_rate].values,
                       "data_group": ["precipitation"]*len(infil[variablesdefinition.infiltration_rate].values),
                       "data_categories": ["infiltration"]*len(infil[variablesdefinition.infiltration_rate].values)
                       }))
    data = data.append(pd.DataFrame({"timestamp": streamflow[variablesdefinition.timestamp].values,
                       "values": streamflow[variablesdefinition.streamflow_label].values,
                       "data_group": ["flow"]*len(streamflow[variablesdefinition.streamflow_label].values),
                       "data_categories": ["measured_flow"]*len(streamflow[variablesdefinition.streamflow_label].values)
                       }))
    data = data.append(pd.DataFrame({"timestamp": mes[variablesdefinition.timestamp].values,
                       "values": mes[variablesdefinition.concentration_label].values,
                       "data_group": ["flow"]*len(mes[variablesdefinition.concentration_label].values),
                       "data_categories": ["computed_flow"]*len(mes[variablesdefinition.concentration_label].values)
                       }))
    print(data)
    grid = sns.FacetGrid(data, col="data_group", hue="data_categories", col_wrap=1, height=1.5)
    grid.map(sns.histplot, "values")
    #grid.map_dataframe(sns.relplot, x="timestamp", y="values")
    # Draw a horizontal line to show the starting point
    #grid.map(plt.axhline, y=0, ls=":", c=".5")

    # Draw a line plot to show the trajectory of each random walk
    #grid.map(plt.plot, "step", "position", marker="o")

    # # Adjust the tick positions and labels
    # grid.set(xticks=np.arange(5), yticks=[-3, 3],
    #          xlim=(-.5, 4.5), ylim=(-3.5, 3.5))

    # Adjust the arrangement of the plots
    grid.fig.tight_layout(w_pad=1)
    # figure(2)
    # #%**********************     Tracé de la pluie      ************************
    # subplot(3,1,1)
    # plot (Pluie.time, Pluie.int/PARAM.dt,'-b')#; hold on      ; % Pluie totale
    # plot (Pluie.time,infil/PARAM.dt,'-k')#; hold off ; % Infiltration
    # title (['MHYDAS UH :   Pluie = ' num2str(L_Pluie) ' mm; Infiltration = ' num2str(L_Inf) ' mm; Ruissellement = ' num2str(L_Ruiss) ' mm']);
    # xlabel ('Temps (HH:MM)');
    # ylabel('Pluie. et Inf. m/s');
    # %legend ('Pluie','Infiltration')
    # datetick('x','HH:MM');
    # hold off
    # grid on
    # #%*************   Tracé des hydrogrammes mesurés et calculés   *************
    # subplot (3,1,2)
    # plot(Pluie.time,Q_cal,'--r');hold on; #% Débit calculé en m3/s
    # plot(streamflow.time,streamflow.int,'-b') ; hold off;  # % Débit mesuré en m3/s
    # title (['Vmes =' num2str(Vol_mes) ' m3; Vcal = ' num2str(Vol_cal) ' m3;    Qmax mes = ' num2str(Qmax_mes) ' m3/s;  Qmax cal = ' num2str(Qmax_cal) ' m3/s; Nash = ' num2str(coeff_Nash) ' ']);
    # xlabel ('Temps (HH:MM)');
    # ylabel ('Débit (m3/s)');
    # %legend ('Débit calculé','Débit mesuré',1)
    # datetick('x','HH:MM');
    # grid on
    # hold off
    # #%********     Tracé des sédimentogrammes mesurés et calculés  *************
    # subplot (3,1,3)
    # plot(Pluie.time,CALC_CONC_TR_LISEM(:,PARAM.Nb_unit+1) ,'--r');hold on; % MES calculé en g/L
    # plot(MES.time,MES.int,'.') ; hold off;   % MES mesuré en g/L
    # title ([' Erosion mesurée = ' num2str(SED_mes) ' kg; Erosion simulée = ' num2str(CALC_Sortie_MES_Parcelle) ' kg; Cmax mes = ' num2str(Cmax_mes) ' g/L;  Cmax cal = ' num2str(Cmax_cal) ' g/L' ]);
    # xlabel ('Temps (HH:MM)');
    # ylabel ('MES (g/L)');
    # %legend ('MES calculé','MES mesuré');
    # datetick('x','HH:MM');
    # #grid on
    # #hold off
    #
    # figure(3)
    # #%*********************     Tableau bilan EROSION     **********************
    # S1 = subplot (1,1,1);
    # title('Bilan sur la parcelle')
    # set(S1, 'XTick', []);
    # set(S1, 'YTick', []);
    # box on
    # Text1{1}= [' BILAN PAR TRONCON ' '**' num2str(metod) ];
    # Text1{end+1}= ['    Production interne  (kg) (' num2str(PARAM.Nb_unit) ' tronçons) : ' ] ;
    # Text1{end+1}= ['    ' num2str(round(CALC_Prod_interne_Tr)) '    (kg)'];
    # Text1{end+1}= ['    ' num2str((CALC_Prod_interne_Tr/(PARAM.Dens_Sed*1000*PARAM.Long_UH*PARAM.Larg_rill*PARAM.Nb_motifs))*1000) '     (mm)'];
    # Text1{end+1}= ' --------------------------------------------';
    # Text1{end+1}= ' BILAN A LA PARCELLE';
    # Text1{end+1}= ['    Erosion mesurée = ' num2str(round(SED_mes)) ' kg <=> ' num2str(10*SED_mes/(PARAM.Long_UH*PARAM.Larg_UH)) ' t/ha' ] ;
    # Text1{end+1}= ['    Erosion simulée =  ' num2str(round(CALC_Sortie_MES_Parcelle)) ' kg <=> ' num2str(10*CALC_Sortie_MES_Parcelle/(PARAM.Long_UH*PARAM.Larg_UH)) ' t/ha'];
    # Text1{end+1}= '                       ';
    # Text1{end+1}= ['    Bilan Erosion diffuse = ' num2str(round(CALC_Splash_Effectif_Parcelle)) ' kg <=> ' num2str(CALC_Splash_Effectif_Parcelle/(PARAM.Dens_Sed*1000*PARAM.Long_UH*PARAM.Larg_UH)*1000) ' mm'] ;
    # Text1{end+1}= ['    Bilan Erosion concentrée = ' num2str(round(sum(CALC_Prod_interne_Tr))) 'kg' ] ;
    # Text1{end+1}= '                       ';
    # Text1{end+1}= ['    Bilan Masse = ' num2str(round(CALC_Sortie_MES_Parcelle-(CALC_Splash_Effectif_Parcelle + sum(CALC_Prod_interne_Tr))) ) 'kg <=> '  num2str(((CALC_Sortie_MES_Parcelle-(CALC_Splash_Effectif_Parcelle + sum(CALC_Prod_interne_Tr)))/CALC_Sortie_MES_Parcelle)*100) ' % de Erosion simulée'];
    # Text1{end+1}= [' Splash total sur sol nu = ' num2str(CALC_Splash_Direct_Tot_Parcelle) ' kg'] ;
    # Text1{end+1}= [' Splash total sous couvert végétal = ' num2str(CALC_Splash_Indirect_Tot_Parcelle) ' kg'] ;
    # Text1{end+1}= [' Coef de Nash = ' num2str(coeff_Nash) ] ;
    #
    # HT1 = text(0.05,0.5,Text1);
    # set(HT1, 'FontSize', 12);
    # set(HT1, 'FontName', 'times new roman');