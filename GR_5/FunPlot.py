import PySimpleGUI as sg
from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GT_arguments as GT_arg;
import GTcomb_arguments as GTcomb_arg
import combustionGT as comb;
import GT2 as gt
import parametricGraphe as pg

"""
The use of this code requires the installation of PySimpleGUI.

This code is used to display windows to allow the user to display the graphs interactively.

WARNING : A LOT OF POINTS MEAN A LOT OF TIME : BE AWARE BEFORE PLOTTING.
"""

# Define the window's contents
layout = [[sg.Text("How many points do you want in the curves?"),sg.Input(key='-INPUT-'),sg.Text(size=(40,1), key='-OUTPUT-')],
          [sg.Text("What graphe do you want to plot?")],
          [sg.Button('Cycle energetic efficiency on the compression ratio with different temperatures')],
          [sg.Button('Cycle energetic efficiency on the cycle work with different temperatures')],
          [sg.Button('Cycle energetic efficiency on the compression ratio with different polytropic efficiency')],
          [sg.Button('Total and cycle energetic efficiency on the compression ratio')],
          [sg.Button('Mechanical work of the cycle on the compression ratio with different temperatures')],
          [sg.Button('Mechanical efficiency on the compression ratio with different temperatures')],
          [sg.Button('All')],
          [sg.Button('Quit')]]

# Create the window
window = sg.Window('Paramatric Graphes', layout)
window2_active = False

# Display and interact with the Window using an Event Loop
while True:
    event, values = window.read()
    # See if user wants to quit or window was closed
    if event == sg.WINDOW_CLOSED or event == 'Quit':
        break
    if event == 'Cycle energetic efficiency on the compression ratio with different temperatures':
        #Check if the user has put an integer as the number of points in the curves
        try:
            val = int(values['-INPUT-'])
            #open a new window to make sure the user wants to plot the graphs even if the code runs for a long time
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                #close the previous window and calls parametricGraphes for plotting the graphics and displaying the progress bar
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('Eta_cyclen_vs_r',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
    if event == 'Cycle energetic efficiency on the compression ratio with different polytropic efficiency':
        try:
            val = int(values['-INPUT-'])
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('Eta_cyclen_vs_r_eta_pic_pit',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
    if event == 'Mechanical work of the cycle on the compression ratio with different temperatures':
        try:
            val = int(values['-INPUT-'])
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('Wcycle_vs_r',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
    if event == 'Mechanical efficiency on the compression ratio with different temperatures':
        try:
            val = int(values['-INPUT-'])
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('eta_mec_vs_r',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
    if event == 'Total and cycle energetic efficiency on the compression ratio':
        try:
            val = int(values['-INPUT-'])
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('eta_cyclen_eta_toten_vs_r',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
    if event == 'Cycle energetic efficiency on the cycle work with different temperatures':
        try:
            val = int(values['-INPUT-'])
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('eta_cyclen_vs_wcy',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
    if event == 'All':
        try:
            val = int(values['-INPUT-'])
            if not window2_active :
                window.Hide()
                window2_active = True
                layout2 = [[sg.Text("The programme may runs for a long time according to the number of points you have chosen")],
                          [sg.Button('Plot'),sg.Button('Exit')]]

                window2 = sg.Window('Paramatric Graphes', layout2)
                while True:
                    event2, values2 = window2.read()
                    if event2 == sg.WINDOW_CLOSED or event2 == 'Exit':
                        window2_active = False
                        window2.Close()
                        window.UnHide()
                        break
                    if event2 == 'Plot':
                        pg.parametricGraphic('all',val)

        except ValueError:
            try:
                val = float(values['-INPUT-'])
                window['-OUTPUT-'].update("The number of points must be an integer")
            except ValueError:
                window['-OUTPUT-'].update("Please indicate an integer")
# Finish up by removing from the screen
window.close()
