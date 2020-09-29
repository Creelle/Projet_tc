
from thermochem import janaf
import matplotlib.pyplot as plt
import numpy as np;

import GT_arguments as GT_arg;

db = janaf.Janafdb();
Mm_CH4 = 0.016; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044
O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
CH4=db.getphasedata('CH4',phase ='g');

def cp_Iconstants(M,T_0,T_1):
    # donne la valeur de cp en J/mol.
    #cp = A + BT + CT^2 + DT^3 (cp évolue en T^3 (sauf au début => à débattre dans le rapport : faire un calcul d'erreur))
    n=1000 # n est arbitraire : j'ai pris n points pour avoir le plus de précision possible
    #db.getphasedata est différente pour l'oxygène que pour les autres gazs donc je dois mettre un if.
    if M=='O2' or M=='N2':
        chim = db.getphasedata(M,'g');
        T = np.linspace(T_0,T_1,n) #je divise en espaces réguliers l intervalle de températures
        cp_values = np.zeros(n)
        for i in range (0,n):
            cp_values[i] = chim.cp(T[i])  #je prends toutes les valeurs de cp correspondant aux températures
        heat_const = np.polyfit(T,cp_values,3) #je trouve les constantes A,B,C,D de de l'interpolation polynomiale grâce à polyfit
        cp_mean =  sum(cp_values)/len(cp_values)
        I = heat_const[3]*(T_1-T_0) + (1/2)*heat_const[2]*(T_1**2-T_0**2) + (1/3)*heat_const[1]*(T_1**3-T_0**3) + (1/4)*heat_const[0]*(T_1**4-T_0**4)
        #j'intègre cp sur l'intervalle de température
    else:
        chim = db.getphasedata(M,phase ='g');
        T = np.linspace(T_0,T_1,n)
        cp_values = np.zeros(n)
        for i in range (0,n):
            cp_values[i] = chim.cp(T[i])
        heat_const = np.polyfit(T,cp_values,3)
        cp_mean =  sum(cp_values)/len(cp_values)
        I = heat_const[3]*(T_1-T_0) + (1/2)*heat_const[2]*(T_1**2-T_0**2) + (1/3)*heat_const[1]*(T_1**3-T_0**3) + (1/4)*heat_const[0]*(T_1**4-T_0**4)
        #je retourne l'intégrale
    return I,cp_mean
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
    #  OUTPUTS
    #  R_f   : ideal gas constant (R*) for specific gas (R/Mm_f) [kJ/kg/K]
    #  m_O2f : mass fraction of oxygen in exhaust gases [-]
    #  m_N2f : mass fraction of nitrogen in exhaust gases [-]
    #  m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
    #  m_H2Of : mass fraction of water steam in exhaust gases [-]
    #  T_out  : outlet gas temperature [K]

    lambda_comb = comb_input.lambda_comb
    x_O2a = comb_input.x_O2a
    x_N2a = comb_input.x_N2a
    coeff = x_N2a/x_O2a
    T_in  = comb_input.T_in +273.15 #[K]
    h_in  = comb_input.h_in #kJ/kg_air
    LHV   =comb_input.LHV #kJ/kg_ch4]
    T_ref =  298.15 #[K]

    # calcul des cp : cec doit encore changer lorsqu on va integrer
    cp_O2 = O2.cp(T_in)
    cp_CO2 = CO2.cp(T_in)
    cp_H2O = H2O.cp(T_in)
    cp_N2 = N2.cp(T_in)
    cp_CH4 = CH4.cp(T_in)



    # CH4 + 2 *lambda * (O2 + coeff*N2) <=> CO2+2*H2O+ 2*lambda*coeff*N2 + 2*(lambda-1)*O2


    Mm_air = x_O2a * Mm_O2 + x_N2a * Mm_N2 # [kg/mol_air]
    ma1 =  Mm_air/Mm_CH4 * 2*lambda_comb/x_O2a # kg_air/kg_CH4 = proportion d air entrant vs combustible
    print('ma1',ma1)
    # T_out =  1400
    # m_a1 = 5
    # iter = 1
    #
    # while iter <= 1000 or error > 0.1 :
    #     Num =
    #     Den1 =
    #     Den2 = (lambda_comb*m_a1 + 1)
    #     T_out_final = 1 + (Num/Den) - h_f0
    #     iter = iter + 1
    #     error = T_out_final - T_out
    #     T_out = T_out_final
    #     print(T_out)
    # print(iter,T_out)


    # A la sortie :
    coeff_stochio = np.array([2*lambda_comb*coeff,1,2,2*(lambda_comb-1)]) # N2- CO2 - H2O - O2
    total_n = sum(coeff_stochio) # nombre de moles total
    molar_conc = coeff_stochio/total_n # concentration des elements mol_co2/mol_t
    molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
    Mm_af = sum(molar_conc*molar_mass) #somme ponderé des masse molaire
    mass_conc = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot

    Cps_out = np.array([cp_N2,cp_CO2,cp_H2O,cp_O2]) # chaleur specifique des gaz sortant J/mol
    H_f0_list = np.array([N2.DeltaH(T_ref),CO2.DeltaH(T_ref),H2O.DeltaH(T_ref),O2.DeltaH(T_ref)])
    # print(H_f0_list)
    hc=(1/Mm_CH4)*cp_Iconstants('CH4',T_ref,T_in)[0] # entalpie sensible du combustibleJ/kg_CH4

    cp_f = np.dot(Cps_out,mass_conc/molar_mass) # valeur du cp des gaz sortant [J/kg/K]

    h_f0 = np.dot(H_f0_list,mass_conc/molar_mass) # valeur d enthalpie de reference des gas sortants [J/kg]

    # print(h_f0)
    T_out = (1 + ((1000*LHV + hc + lambda_comb*ma1*h_in*1000)/((lambda_comb*ma1+1)*cp_f*T_ref)) - h_f0/(cp_f*T_ref))*T_ref
    # E_left = (LHV*1000+hc+lambda_comb*ma1*h_in*1000)*(T_in-T_ref)
    # E_right = Cp_f*ma1
    # T_out = E_left/(E_right*30)
    outputs = GT_arg.comb_output();
    outputs.T_out = T_out
    outputs.R_f = 8.31/1000/Mm_af # [kJ/kg/K]
    outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = mass_conc  #[-]
    return outputs;

sol =combustionGT(GT_arg.comb_input(lambda_comb = 2))
print(sol.T_out)
