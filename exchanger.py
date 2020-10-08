from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GT_arguments as GT_arg;
import combustionGT as comb;

O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
Mm_O2 = 0.032;#kg/mol
Mm_N2 = 0.028;#kg/mol
conc_O2 = 0.21;# 21% in molar
conc_N2 = 0.79;# 79% in molar


def heatexchanger(exchanger_input,T_out):
    """
    On connait le mass flow d'air in et le mass flow de fumée in ainsi que tous les cp correspondants
    --> on les connait car on veut un certains amount of power à la sortie directement lié au massflow dans la GT.
    On connait la température in de l'air et la température in des fumées (calculée par combustionGT).
    J'indique dans la défintion de la fonction, la température out de l'air que je souhaite à la sortie de l'échangeur
    Return la température des fumées à la sortie de l'échangeur ainsi que la surface d'échange nécessaire (donné grâce à Hausbrand)
    """

    U = exchanger_input.U
    Mflow_air_in = exchanger_input.Mflow_air_in
    Mflow_f_in = exchanger_input.Mflow_f_in
    T_air_in = exchanger_input.T_air_in
    T_f_in = exchanegr_input.T_f_in







    return outputs
