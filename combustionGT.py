
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

def cp_air(T):#kJ/kg/K
    Mm_a = 0.21 * 0.032 + 0.79 * 0.028;
    m_O2 = (0.21*0.032)/Mm_a;# mass proportion of O2
    m_N2 = (0.79*0.028)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg_air/K
    #R = 8.31/Mm_a/1000
    #gamma = Cp/(Cp-R)
    return Cp;
def cp_air_T(T):#kJ/kg/K
    Mm_a = 0.21 * 0.032 + 0.79 * 0.028;
    m_O2 = (0.21*0.032)/Mm_a;# mass proportion of O2
    m_N2 = (0.79*0.028)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg_air/K
    #R = 8.31/Mm_a/1000
    #gamma = Cp/(Cp-R)
    return Cp/T;

def cpCH4(T):
    return CH4.cp(T)
def cpCH4_T(T):
    return CH4.cp(T)*(1/T)
def cpN2(T):
    return N2.cp(T)
def cpO2(T):
    return O2.cp(T)
def cpH2O(T):
    return H2O.cp(T)
def cpCO2(T):
    return CO2.cp(T)
def janaf_integrate(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt) # int(cp)dt [J/mol]]
def cp_mean(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)/len(values)) #  cp_mean [J/mol/K]
# T1= 20
# T2 = 500
# print('janaf_integrate',janaf_integrate(cpN2,T1,T2,0.0001))
# print('cpIconstants',cp_Iconstants('N2',T1,T2))

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
    HHV = comb_input.HHV #kJ/kg_CH4
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

    h_f0 = np.dot([N2.hef(T0),CO2.hef(T0),H2O.hef(T0),O2.hef(T0)],mass_conc/molar_mass)*1000 # valeur d enthalpie sensible des gas sortants

    # hc=(1/Mm_CH4)*cp_Iconstants('CH4',T0-15,T_in)[0] # entalpie sensible du combustibleJ/kg_CH4
    # hc2=(1/Mm_CH4)*(CH4.hef(T_in)-CH4.hef(T0-15))*1000 # 2 maniere
    hc= janaf_integrate(cpCH4,T0-15,T_in,0.0001)/Mm_CH4 # 3 maniere
    # print('hello', hc,hc2,hc3)

    #  calcul de T_out par iteration
    iter = 1
    T_out = 2000 #valeur initiale
    error = 1
    dt = 0.1

    while iter < 1000 and error > 0.01 :
        #cps_out = np.array([cp_Iconstants('N2',T0,T_out)[1],cp_Iconstants('CO2',T0,T_out)[1],cp_Iconstants('H2O',T0,T_out)[1],cp_Iconstants('O2',T0,T_out)[1]])
        cps_out = np.array([cp_mean(cpN2,T0,T_out,dt),cp_mean(cpCO2,T0,T_out,dt),cp_mean(cpH2O,T0,T_out,dt),cp_mean(cpO2,T0,T_out,dt)])
        cp_f = np.dot(cps_out,mass_conc/molar_mass) #J/kg/K
        T_out_final = (T0 + ((1000*LHV + hc + lambda_comb*ma1*h_in*1000)/((lambda_comb*ma1+1)*cp_f)) - h_f0/(cp_f))
        iter = iter + 1
        error = abs(T_out_final - T_out)
        T_out = T_out_final
    print("Nombre d'itérations : ",iter)
    print("T_out : ",T_out,"K")

    #calcul de l exergie et eta_combex ==> see formula page 28
    e_c = HHV+ 15 * (cp_mean(cpCH4,273.15,T0,dt)/0.016 + mass_conc[3]/0.032*cp_mean(cpO2,273.15,T0,dt) - mass_conc[1]/0.044*cp_mean(cpCO2,273.15,T0,dt) - mass_conc[2]/0.018*cp_mean(cpH2O,273.15,T0,dt))/1000 - T0*(CH4.S(273.15)/0.016+cp_mean(cpCH4,273.15,T0,dt)/0.016*np.log(T0/273.15))/1000 # kJ/kg_ch4
    #e_c2 = HHV-T0*CH4.S(273.15)/0.016/1000 voir page 28
    hf = (LHV*1000+ hc + lambda_comb*ma1*h_in*1000)/( lambda_comb*ma1+1) # a verifier J/kg

    # H_f_list = np.array([N2.hef(T_out),CO2.hef(T_out),H2O.hef(T_out),O2.hef(T_out)])*1000
    # h_f2 = np.dot(H_f_list,mass_conc/molar_mass)
    Sf = np.dot([N2.S(T_out),CO2.S(T_out),H2O.S(T_out),O2.S(T_out)],mass_conc/molar_mass) #J/kg/K N2- CO2 - H2O - O2
    Sf0 = np.dot([N2.S(T0),CO2.S(T0),H2O.S(T0),O2.S(T0)],mass_conc/molar_mass)
    e_f= (hf-h_f0 -T0*(Sf-Sf0))/1000#exergie des fumées kJ/kg_CH4

    """
    Version 1 de e_a

    heat_const_O2 = cp_Iconstants('O2',T0,T_in)[2]
    heat_const_N2 = cp_Iconstants('N2',T0,T_in)[2]
    cp_bar_O2 = heat_const_O2[3]*np.log(T_in/T0) + heat_const_O2[2]*(T_in-T0) + (1/2)*heat_const_O2[1]*(T_in**2-T0**2) + (1/3)*heat_const_O2[0]*(T_in**3-T0**3)
    cp_bar_N2 = heat_const_N2[3]*np.log(T_in/T0) + heat_const_N2[2]*(T_in-T0) + (1/2)*heat_const_N2[1]*(T_in**2-T0**2) + (1/3)*heat_const_N2[0]*(T_in**3-T0**3)
    e_a_2 = ((0.21/0.032*cp_mean(cpO2,T0,T_in,dt) + 0.79/0.028*cp_mean(cpN2,T0,T_in,dt))*(T_in-T0) - (0.21/0.032*cp_bar_O2 + 0.79/0.028*cp_bar_N2)*T0)/1000
    print("e_a_2 =", e_a_2)
    """
    e_a = cp_mean(cp_air,T0,T_in,dt)*(T_in-T0) - janaf_integrate(cp_air_T,T0,T_in,dt)*T0
    print("e_a =", e_a)
    e_cr = ((cp_mean(cpCH4,T0,T_in,dt)/0.016)*(T_in-T0) - janaf_integrate(cpCH4_T,T0,T_in,dt)*T0/0.016)/1000 #kJ/kg_CH4
    e_r = e_a*((lambda_comb*ma1)/(lambda_comb*ma1+1)) + e_cr*(1/(lambda_comb*ma1+1))
    eta_combex = (e_f-e_r)*(lambda_comb*ma1+1)/e_c# voir page 31 #bon il y a un probleme obviously
    print("eta_combex = ",eta_combex)

    # remplissage des outputs
    outputs = GT_arg.comb_output();
    outputs.T_out = T_out
    outputs.R_f = 8.31/1000/Mm_af # [kJ/kg/K]
    outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = mass_conc  #[-]
    outputs.Mm_af = Mm_af
    outputs.lambda_comb = lambda_comb
    outputs.ma1 = ma1
    return outputs;

sol =combustionGT(GT_arg.comb_input(lambda_comb = 5))
#print(sol.T_out)

"""
#Fais le plot de T_out en fonction de lambda_comb

x = np.linspace(1,10,5)
y = np.zeros(len(x))
for i in range (0,len(x)) :
    y[i] = combustionGT(GT_arg.comb_input(lambda_comb = x[i])).T_out
plt.plot(x,y)
plt.ylabel('Temperature [K]')
plt.xlabel('Lambda')
plt.show()
"""

#Mettre
#combustionGT --> return eta_combex
# x = np.linspace(1,10,10)
# y = np.zeros(len(x))
# for i in range (0,len(x)) :
#     y[i] = combustionGT(GT_arg.comb_input(lambda_comb = x[i]))
# plt.plot(x,y)
# plt.ylabel('etacombex [K]')
# plt.xlabel('Lambda')
# plt.show()
