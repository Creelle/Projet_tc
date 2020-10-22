from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GT_arguments as GT_arg;
import GTcomb_arguments as GTcomb_arg
import combustionGT as comb;
import GT as gt

"""
J'ai crée une fonction pour afficher les graphes pour le moment 3 choix possibles :
1)Le graphique de eta cyclen sur le taux de compression au compresseur = 'Eta_cyclen_vs_r'
2)Le graphique du travail moteur sur le taux de compression au compresseur = 'Wcycle_vs_r'
3)Tous les graphiques = 'all'

Les arguments de la fonction parametricGraphic sont donc M = {'Eta_cyclen_vs_r','Wcycle_vs_r','all'} et number = le nombre de points par courbe (il y a 4 courbes à chaque fois)

ATTENTION : estimation de temps de run de la fonction =~ 1sec par point!!!!!


"""

def parametricGraphic(M,number):
    if M == 'Eta_cyclen_vs_r' or M == 'all':
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        for j in range (0,4):
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = 1273.15)).eta[0]
            if j==1:
                for i in range (0,len(x)) :
                    y2[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = 1473.15)).eta[0]
            if j==2:
                for i in range (0,len(x)) :
                    y3[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = 1673.15)).eta[0]
            if j==3:
                for i in range (0,len(x)) :
                    y4[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = 1873.15)).eta[0]

        fig=plt.figure()
        plt.plot(x,y1,x,y2,x,y3,x,y4)
        plt.legend(('T = 1000 °C', 'T = 1200 °C', 'T = 1400 °C','T = 1600 °C'),
                   loc='upper right')
        plt.ylabel('Eta cyclen')
        plt.xlabel('r')
        plt.savefig('figures/Eta_cyclen_vs_r.png')

    if M == 'Wcycle_vs_r' or M == 'all':
        T=np.array([1073.15,1273.15,1473.15,1673.15,1873.15])
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        y5 = np.zeros(len(x))
        for j in range (0,5):
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
            if j==1:
                for i in range (0,len(x)) :
                    y2[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
            if j==2:
                for i in range (0,len(x)) :
                    y3[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
            if j==3:
                for i in range (0,len(x)) :
                    y4[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
            if j==4:
                for i in range (0,len(x)) :
                    y5[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
        plt.plot(x,y1,x,y2,x,y3,x,y4,x,y5)
        plt.legend(('T = 800 °C','T = 1000 °C', 'T = 1200 °C', 'T = 1400 °C','T = 1600 °C'),
                   loc='upper right')
        plt.ylabel('Wmcy')
        plt.xlabel('r')
        plt.show()

    """
    x = np.linspace(1.1,100,20)
    y = np.zeros(len(x))
    for i in range (0,len(x)) :
        y[i] = gt.GT_simple(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=288.15,T_0 = 288.15,r=x[i],k_cc=0.95,T3 = 1273.15)).eta[0]
    plt.plot(x,y)
    plt.ylabel('Eta cyclen')
    plt.xlabel('r')
    plt.show()
    """
parametricGraphic('Eta_cyclen_vs_r',2)
