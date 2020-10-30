from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;
import PySimpleGUI as sg

import GT_arguments as GT_arg;
import GTcomb_arguments as GTcomb_arg
import combustionGT as comb;
import GT2 as gt
plt.rcParams.update({'font.size': 10})
"""
J'ai crée une fonction pour afficher les graphes pour le moment 3 choix possibles :
1)Le graphique de eta cyclen sur le taux de compression au compresseur = 'Eta_cyclen_vs_r'
2)Le graphique du travail moteur sur le taux de compression au compresseur = 'Wcycle_vs_r'
3)Tous les graphiques = 'all'

Les arguments de la fonction parametricGraphic sont donc M = {'Eta_cyclen_vs_r','Eta_cyclen_vs_r_eta_pic_pit','Wcycle_vs_r','eta_mec_vs_r','eta_cyclen_eta_toten_vs_r','eta_cyclen_vs_wcy','all'} et number = le nombre de points par courbe (il y a 4 courbes à chaque fois)

ATTENTION : estimation de temps de run de la fonction =~ 1sec par point!!!!!


"""

def parametricGraphic(M,number):
    if M == 'Eta_cyclen_vs_r' or M == 'all':
        T=np.array([1000,1200,1400,1600])#[°C]
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        count = 0
        # layout the form
        layout3 = [[sg.Text('\u03B7 cyclen according to r and temperatures')],
                  [sg.ProgressBar(1, orientation='h', size=(20,20), key='progress')],
                  [sg.Cancel()]]

        # create the form
        window3 = sg.Window('Paramatric Graphes Progress Meter').Layout(layout3)
        progress_bar = window3.FindElement('progress')

        for j in range (0,4):
            # check to see if the cancel button was clicked and exit loop if clicked
            event3, values3 = window3.Read(timeout=0)
            if event3 == 'Cancel' or event3 == None:
                break
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==1:
                for i in range (0,len(x)) :
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
                    y2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==2:
                for i in range (0,len(x)) :
                    y3[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==3:
                for i in range (0,len(x)) :
                    y4[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)

        fig0=plt.figure()
        plt.plot(x,y1,x,y2,x,y3,x,y4)
        plt.legend(('t3 = 1000 °C', 't3 = 1200 °C', 't3 = 1400 °C','t3 = 1600 °C'),
                   loc='upper left')
        plt.title('Cycle energetic efficiency on the compression ratio \n with different temperatures t3. \n Pe = 230 MW, t1 = 15 °C, k_mec = 0.015, k_cc = 0.95, \u03B7pC = \u03B7pT = 0.90')
        plt.ylabel('\u03B7 cyclen')
        plt.xlabel('r')
        plt.savefig('figures/Eta_cyclen_vs_r.png')
        window3.Close()

    if M == 'Eta_cyclen_vs_r_eta_pic_pit' or M == 'all':
        eta_pic_pit=np.array([0.86,0.88,0.9,0.92])#[°C]
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        count = 0
        # layout the form
        layout3 = [[sg.Text('\u03B7 cyclen according to r and \u03B7 PiC = \u03B7 PiT')],
                  [sg.ProgressBar(1, orientation='h', size=(20,20), key='progress')],
                  [sg.Cancel()]]

        # create the form
        window3 = sg.Window('Paramatric Graphes Progress Meter').Layout(layout3)
        progress_bar = window3.FindElement('progress')
        for j in range (0,4):
            # check to see if the cancel button was clicked and exit loop if clicked
            event3, values3 = window3.Read(timeout=0)
            if event3 == 'Cancel' or event3 == None:
                break
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = 1400,eta_PiC=eta_pic_pit[j],eta_PiT=eta_pic_pit[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==1:
                for i in range (0,len(x)) :
                    y2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = 1400,eta_PiC=eta_pic_pit[j],eta_PiT=eta_pic_pit[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==2:
                for i in range (0,len(x)) :
                    y3[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = 1400,eta_PiC=eta_pic_pit[j],eta_PiT=eta_pic_pit[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==3:
                for i in range (0,len(x)) :
                    y4[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = 1400,eta_PiC=eta_pic_pit[j],eta_PiT=eta_pic_pit[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)

        fig1=plt.figure()
        plt.plot(x,y1,x,y2,x,y3,x,y4)
        plt.legend(('\u03B7pC = \u03B7pT = 0.86', '\u03B7pC = \u03B7pT = 0.88', '\u03B7pC = \u03B7pT = 0.9','\u03B7pC = \u03B7pT = 0.92'),
                   loc='upper left')
        plt.ylabel('\u03B7 cyclen')
        plt.xlabel('r')
        plt.title('Cycle energetic efficiency on the compression ratio \n with different polytropic efficiency \n Pe = 230 MW, t1 = 15 °C, t3 = 1400°C, k_mec = 0.015, k_cc = 0.95')
        plt.savefig('figures/Eta_cyclen_vs_r_eta_pic_pit.png')
        window3.Close()

    if M == 'Wcycle_vs_r' or M == 'all':
        T=np.array([1000,1200,1400,1600])#[°C]
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        count = 0
        # layout the form
        layout3 = [[sg.Text('Wcy according to r and temperatures')],
                  [sg.ProgressBar(1, orientation='h', size=(20,20), key='progress')],
                  [sg.Cancel()]]

        # create the form
        window3 = sg.Window('Paramatric Graphes Progress Meter').Layout(layout3)
        progress_bar = window3.FindElement('progress')
        for j in range (0,4):
            # check to see if the cancel button was clicked and exit loop if clicked
            event3, values3 = window3.Read(timeout=0)
            if event3 == 'Cancel' or event3 == None:
                break
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==1:
                for i in range (0,len(x)) :
                    y2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==2:
                for i in range (0,len(x)) :
                    y3[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==3:
                for i in range (0,len(x)) :
                    y4[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
        fig2=plt.figure()
        plt.plot(x,y1,x,y2,x,y3,x,y4)
        plt.legend(('t3 = 1000 °C', 't3 = 1200 °C', 't3 = 1400 °C','t3 = 1600 °C'),
                   loc='upper left')
        plt.ylabel('Wmcy [kJ/kg]')
        plt.xlabel('r')
        plt.title('Mechanical work of the cycle on the compression ratio \n with different temperatures t3 \n Pe = 230 MW, t1 = 15 °C, k_mec = 0.015, k_cc = 0.95, \u03B7pC = \u03B7pT = 0.90')
        plt.savefig('figures/Wcycle_vs_r.png')
        window3.Close()
    if M == 'eta_mec_vs_r' or M == 'all':
        T=np.array([1000,1200,1400,1600])#[°C]
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        count = 0
        # layout the form
        layout3 = [[sg.Text('\u03B7 mec according to r and temperatures')],
                  [sg.ProgressBar(1, orientation='h', size=(20,20), key='progress')],
                  [sg.Cancel()]]

        # create the form
        window3 = sg.Window('Paramatric Graphes Progress Meter').Layout(layout3)
        progress_bar = window3.FindElement('progress')
        for j in range (0,4):
            # check to see if the cancel button was clicked and exit loop if clicked
            event3, values3 = window3.Read(timeout=0)
            if event3 == 'Cancel' or event3 == None:
                break
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[2]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==1:
                for i in range (0,len(x)) :
                    y2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[2]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==2:
                for i in range (0,len(x)) :
                    y3[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[2]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
            if j==3:
                for i in range (0,len(x)) :
                    y4[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[2]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*4)
        fig3=plt.figure()
        plt.plot(x,y1,x,y2,x,y3,x,y4)
        plt.legend(('t3 = 1000 °C', 't3 = 1200 °C', 't3 = 1400 °C','t3 = 1600 °C'),
                   loc='upper left')
        plt.ylabel('\u03B7mec')
        plt.xlabel('r')
        plt.title('Mechanical efficiency on the compression ratio \n with different temperatures t3 \n Pe = 230 MW, t1 = 15 °C, k_mec = 0.015, k_cc = 0.95, \u03B7pC = \u03B7pT = 0.90')
        plt.savefig('figures/eta_mec_vs_r.png')
        window3.Close()
    if M == 'eta_cyclen_eta_toten_vs_r' or M == 'all':
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        count = 0
        # layout the form
        layout3 = [[sg.Text('\u03B7 cyclen and \u03B7 toten according to r')],
                  [sg.ProgressBar(1, orientation='h', size=(20,20), key='progress')],
                  [sg.Cancel()]]

        # create the form
        window3 = sg.Window('Paramatric Graphes Progress Meter').Layout(layout3)
        progress_bar = window3.FindElement('progress')
        for j in range (0,2):
            # check to see if the cancel button was clicked and exit loop if clicked
            event3, values3 = window3.Read(timeout=0)
            if event3 == 'Cancel' or event3 == None:
                break
            if j==0:
                for i in range (0,len(x)) :
                    y1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = 1400))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*2)
            if j==1:
                for i in range (0,len(x)) :
                    y2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = 1400))[3]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*2)
        fig4=plt.figure()
        plt.plot(x,y1,x,y2)
        plt.legend(('\u03B7 cyclen', '\u03B7 toten'),
                   loc='upper left')
        plt.ylabel('\u03B7')
        plt.xlabel('r')
        plt.title('Total and cycle energetic efficiency on the compression ratio \n Pe = 230 MW, t1 = 15 °C, t3 = 1400°C, k_mec = 0.015, \n k_cc = 0.95, \u03B7pC = \u03B7pT = 0.90')
        plt.savefig('figures/eta_cyclen_eta_toten_vs_r.png')
        window3.Close()
    if M == 'eta_cyclen_vs_wcy' or M == 'all':
        T=np.array([1000,1200,1400,1600])#[°C]
        x = np.linspace(2,100,number)
        y1 = np.zeros(len(x))
        y2 = np.zeros(len(x))
        y3 = np.zeros(len(x))
        y4 = np.zeros(len(x))
        x1 = np.zeros(len(x))
        x2 = np.zeros(len(x))
        x3 = np.zeros(len(x))
        x4 = np.zeros(len(x))
        count = 0
        # layout the form
        layout3 = [[sg.Text('\u03B7 cyclen according to Wcy and temperatures')],
                  [sg.ProgressBar(1, orientation='h', size=(20,20), key='progress')],
                  [sg.Cancel()]]

        # create the form
        window3 = sg.Window('Paramatric Graphes Progress Meter').Layout(layout3)
        progress_bar = window3.FindElement('progress')
        for j in range (0,4):
            # check to see if the cancel button was clicked and exit loop if clicked
            event3, values3 = window3.Read(timeout=0)
            if event3 == 'Cancel' or event3 == None:
                break
            if j==0:
                for i in range (0,len(x)) :
                    x1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    y1[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*8)
            if j==1:
                for i in range (0,len(x)) :
                    x2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    y2[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*8)
            if j==2:
                for i in range (0,len(x)) :
                    x3[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    y3[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*8)
            if j==3:
                for i in range (0,len(x)) :
                    x4[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[1]
                    y4[i] = gt.GT(GT_arg.GT_input(Pe = 230e3,k_mec =0.015, T_ext=15,T_0 = 15,r=x[i],k_cc=0.95,T3 = T[j]))[0]
                    count =count +1
                    progress_bar.UpdateBar(count+1, number*8)
        fig5=plt.figure()
        plt.plot(x1,y1,x2,y2,x3,y3,x4,y4)
        plt.legend(('t3 = 1000 °C', 't3 = 1200 °C', 't3 = 1400 °C','t3 = 1600 °C'),
                   loc='upper left')
        plt.ylabel('\u03B7 cyclen')
        plt.xlabel('W_cycle [kJ/kg]')
        plt.title('Cycle energetic efficiency on the cycle work \n with different temperatures t3 \n Pe = 230 MW, t1 = 15 °C, k_mec = 0.015, k_cc = 0.95, \u03B7pC = \u03B7pT = 0.90')
        plt.savefig('figures/eta_cyclen_vs_wcy.png')
        window3.Close()
    plt.show()
#parametricGraphic('Eta_cyclen_vs_r_eta_pic_pit',2)
