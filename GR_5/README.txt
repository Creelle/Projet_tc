
Title :  GT project 2020 in LMECA2150 Thermal cycles

Objectives :
  Analyse the performances of a simple gas turbine from an energetic and exergetic
  point of view

Files:
  The main file is GT.py. This file contains the function GT() which computes
  the outputs of a gas turbine from the parameters given in input.
  These parameters are given in the file GT_arguments.py

  The gas turbine has a combustion chamber. The file combustionGT.py computes
  the exit of the combustion chamber  from inputs. It also makes an exergetic analysis.
  The inputs and outputs are defined in the file GTcomb_arguments

  Also in GTcomb_arguments, inputs and outputs for a heat exchanger are defined.
  This heat exchanger file: exchanger.py is only functional for two flows of
  air but a more general version can be considered in latter versions.

  Anyway, the file exchanger.py takes care of those input and output arguments.
  In GTexchanger.py, a heat exchanger is added in the cycle between the compressor
  and the combustion chamber.

  Besides, the files parametricGraphe.py, GT2.py and windowtest.py are there to
  make fancy graphes of the eta_cyclen, eta_toten, eta_mec and Wm as a function
  of several inputs. To be more precise:
  - GT2.py : a lighter version of GT which only returns the outputs needed
  - parametricGraphe.py:  the beast file where the loops are made to compute
                          with the different parameters
  - windowtest.py : this file permits to do a fancy interface with user when he
                    wants to plot a parametric graph from parametricGraphe.py

  Last but not least, useful.py is the backspine of our files as it contains
  all the useful functions used in GT, combustionGT, exchanger,Gtexchanger 
