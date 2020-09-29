
from thermochem import janaf
from thermochem import combustion
from thermochem import burcat
from numpy import*
from scipy.linalg import*
import matplotlib.pyplot as plt
db = janaf.Janafdb();
burcatdb = burcat.Elementdb()
# mix = burcatdb.getmixturedata([("O2 REF ELEMENT", 20.9476), ("N2 REF ELEMENT",78.084), ("CO2", 0.0319), ("AR REF ELEMENT", 0.9365), ])

import numpy as np;

import GT_arguments as GT_arg;

def cp_Iconstants(M,T_0,T_1):
    #cp = A + BT + CT^2 + DT^3 (cp évolue en T^3 (sauf au début => à débattre dans le rapport : faire un calcul d'erreur))
    n=1000 # n est arbitraire : j'ai pris n points pour avoir le plus de précision possible
    #db.getphasedata est différente pour l'oxygène que pour les autres gazs donc je dois mettre un if.
    if M=='O2' or M=='N2':
        chim = db.getphasedata(M,'g');
        T = linspace(T_0,T_1,n) #je divise en espaces réguliers l intervalle de températures
        cp_values = zeros(n)
        for i in range (0,n):
            cp_values[i] = chim.cp(T[i])  #je prends toutes les valeurs de cp correspondant aux températures
        heat_const = polyfit(T,cp_values,3) #je trouve les constantes A,B,C,D de de l'interpolation polynomiale grâce à polyfit
        I = heat_const[3]*(T_1-T_0) + (1/2)*heat_const[2]*(T_1**2-T_0**2) + (1/3)*heat_const[1]*(T_1**3-T_0**3) + (1/4)*heat_const[0]*(T_1**4-T_0**4)
        #j'intègre cp sur l'intervalle de température
    else:
        chim = db.getphasedata(M,phase ='g');
        T = linspace(T_0,T_1,n)
        cp_values = zeros(n)
        for i in range (0,n):
            cp_values[i] = chim.cp(T[i])
        heat_const = polyfit(T,cp_values,3)
        I = heat_const[3]*(T_1-T_0) + (1/2)*heat_const[2]*(T_1**2-T_0**2) + (1/3)*heat_const[1]*(T_1**3-T_0**3) + (1/4)*heat_const[0]*(T_1**4-T_0**4)
        #je retourne l'intégrale
    return I
print("O2 :", cp_Iconstants('O2',273,600))

def combustionGT(comb_input):
    """
     GT Gas turbine modelisation
     GT(P_e,options,display) compute the thermodynamics states for a Gas
     turbine based on several inputs (given in OPTION) and based on a given
     electricity production P_e. It returns the main results. It can as well
     plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)

     INPUTS (some inputs can be dependent on others => only one of these 2 can
             be activated) Refer to Fig 3.1 from reference book (in english)0
"""
    #fuel=burcatdb.getelementdata("CH4   RRHO")  # choose the right fuel
    #combObject=combustion.SimpleCombustor(fuel,1/lambda_comb,burcatdb)

    # donne une liste avec les coeff stochio pour les produits de la combustion dans cette ordre N2 - CO2 - H2O - O2
    #coeff_stochio = combustion.balance(fuel,1,1/lambda_comb)[1]
    #print(combustion.balance(fuel,1,1/lambda_comb))
    #list_coeff = np.asarray(list(coeff_stochio.values()))

    #total_n = sum(list_coeff) # nombre de moles total
    #molar_conc = list_coeff/total_n # concentration des elements mol_co2/mol_t
    #molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
    #Mm_af = sum(molar_conc*molar_mass) #somme ponderé des masse molaire

    #  OUTPUTS
    #  R_f   : ideal gas constant (R*) for specific gas (R/Mm_f) [kJ/kg/K]
    #  m_O2f : mass fraction of oxygen in exhaust gases [-]
    #  m_N2f : mass fraction of nitrogen in exhaust gases [-]
    #  m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
    #  m_H2Of : mass fraction of water steam in exhaust gases [-]
    #  T_out  : outlet gas temperature [K]

    #outputs = GT_arg.comb_output();
    #outputs.T_out = combObject.adiabatic_flame_temp(273.15+T_in); #[K]
    #outputs.R_f = 8.31/1000/Mm_af # [kJ/kg/K]
    #outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = molar_conc*molar_mass/Mm_af #[-]
    Mm_CH4 = 0.016; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044
    T_out =  1400
    m_a1 = 5
    iter = 1

    while iter <= 1000 or error > 0.1 :
        Num =
        Den1 =
        Den2 = (lambda_comb*m_a1 + 1)
        T_out_final = 1 + (Num/Den) - h_f0
        iter = iter + 1
        error = T_out_final - T_out
        T_out = T_out_final
        print(T_out)
    print(iter,T_out)





    outputs = 0;
    #outputs.T_out = T_out
    #outputs.R_f = 8.31/1000 # [kJ/kg/K]
    #outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = 3 #[-]
    return outputs;
### test
#print()
#sol =combustionGT(GT_arg.comb_input)
#print(sol.T_out)
#combustionGT(GT_arg.comb_input(lambda_comb = 2))
