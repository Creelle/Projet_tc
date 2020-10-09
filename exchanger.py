from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GT_arguments as GT_arg;
import combustionGT as comb;

db = janaf.Janafdb();
Mm_CH4 = 0.016; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044
O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
CH4=db.getphasedata('CH4',phase ='g');
CO = db.getphasedata('CO',phase ='g')

def cp_air(T):#kJ/kg/K ---> OK
    Mm_a = 0.21 * 0.032 + 0.79 * 0.028;
    m_O2 = (0.21*0.032)/Mm_a;# mass proportion of O2
    m_N2 = (0.79*0.028)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    cp = cp_a/Mm_a;#J/kg_air/K
    #R = 8.31/Mm_a/1000
    #gamma = Cp/(Cp-R)
    return cp;
def cp_air_T(T):#kJ/kg/K
    Mm_a = 0.21 * 0.032 + 0.79 * 0.028;
    m_O2 = (0.21*0.032)/Mm_a;# mass proportion of O2
    m_N2 = (0.79*0.028)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    cp = cp_a/Mm_a;#J/kg_air/K
    #R = 8.31/Mm_a/1000
    #gamma = Cp/(Cp-R)
    return cp/T;

def cpCH4(T):
    return CH4.cp(T)
def cpCH4_T(T):
    return CH4.cp(T)*(1/T)
def cpO2(T):
    return O2.cp(T)
def cpO2_T(T):
    return O2.cp(T)*(1/T)
def cpH2O(T):
    return H2O.cp(T)
def cpH2O_T(T):
    return H2O.cp(T)*(1/T)
def cpCO2(T):
    return CO2.cp(T)
def cpCO2_T(T):
    return CO2.cp(T)*(1/T)
def cpN2(T):
    return N2.cp(T)
def cpN2_T(T):
    return N2.cp(T)*(1/T)
def janaf_integrate(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt) # int(cp)dt [J/mol/K]]
def cp_mean(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)/len(values)) #  cp_mean [J/mol/K]

U = exchanger_input.U
Mflow_air_in = exchanger_input.Mflow_air_in
Mflow_f_in = exchanger_input.Mflow_f_in
T_air_in = exchanger_input.T_air_in
T_f_in = exchanger_input.T_f_in
lambda_comb = exchanger_input.comb_lambda

x_O2a = 0.21
x_N2a = 0.79
coeff = x_N2a/x_O2a
T0 =  288.15 #[K]

# CH4 + 2 *lambda * (O2 + coeff*N2) <=> CO2+2*H2O+ 2*lambda*coeff*N2 + 2*(lambda-1)*O2



Mm_air = x_O2a * Mm_O2 + x_N2a * Mm_N2 # [kg/mol_air]
ma1 =  Mm_air/Mm_CH4 * 2/x_O2a # kg_air/kg_CH4 = proportion d air entrant vs combustible

# A la sortie :
coeff_stochio = np.array([2*lambda_comb*coeff,1,2,2*(lambda_comb-1)]) # N2- CO2 - H2O - O2
total_n = sum(coeff_stochio) # nombre de moles total
molar_conc = coeff_stochio/total_n # concentration des elements mol_co2/mol_t
molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
Mm_af = sum(molar_conc*molar_mass) #somme ponderé des masse molaire
mass_conc = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot

def heatexchanger(exchanger_input,T_air_out):
    """
    On connait le mass flow d'air in et le mass flow de fumée in ainsi que tous les cp correspondants
    --> on les connait car on veut un certains amount of power à la sortie directement lié au massflow dans la GT.
    On connait la température in de l'air et la température in des fumées (calculée par combustionGT).
    J'indique dans la défintion de la fonction, la température out de l'air que je souhaite à la sortie de l'échangeur
    Return la température des fumées à la sortie de l'échangeur ainsi que la surface d'échange nécessaire (donné grâce à Hausbrand)
    """
    dt = 0.01
    iter = 1
    T_f_out = 1000 #valeur initiale
    print(T_f_out)
    error = 1

    outputs = GT_arg.exchanger_output();
    outputs.T_air_out = T_air_out

    Q = Mflow_air_in * janaf_integrate(cp_air,T_air_in,T_air_out,dt)

    while iter<1000 and errror>0.1 :
        cps_out = np.array([janaf_integrate(cpN2,T0,T_out,dt),janaf_integrate(cpCO2,T0,T_out,dt),janaf_integrate(cpH2O,T0,T_out,dt),janaf_integrate(cpO2,T0,T_out,dt)])
        cp_f = np.dot(cps_out,mass_conc/molar_mass) #J/kg_fumée/K


    outputs.Q = Q

    return outputs

sol = heatexchanger(GT_arg.exchanger_input(U = 2),600)
print(sol.Q)
