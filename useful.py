from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;

"""
fichier contenant les fonctions utiles pour combustionGT et GT
"""

O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
CH4=db.getphasedata('CH4',phase ='g');
Mm_O2 = 0.032;#kg/mol
Mm_N2 = 0.028;#kg/mol
conc_O2 = 0.21;# 21% in molar
conc_N2 = 0.79;# 79% in molar

def air_mixture(T):#kJ/kg/K
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg/K
    R = 8.31/Mm_a/1000
    gamma = Cp/(Cp-R)
    return Cp,gamma,R;

def cp_air(T,conc_mass,Mm_a):
    cps = np.array([N2.cp(T),CO2.cp(T),H2O.cp(T),O2.cp(T)])
    cp_air = np.dot(conc_mass,cps);#J/mol/K
    return cp_air/Mm_a #J/kg/K

def cpCH4(T):
    return CH4.cp(T)
def cpCH4_T(T):
    return CH4.cp(T)*(1/T)
def cpO2(T):
    return O2.cp(T)
def cpCO2(T):
    return CO2.cp(T)
def cpN2(T):
    return N2.cp(T)
def cpH2O(T):
    return H2O.cp(T)
def cp_mean(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)/len(values)) #  cp_mean [J/mol/K]


#fonction qui donne l enthalpie (kJ/kg), T temperature concentration massique mass : array (N2 , CO2, H20,O2)
#molar mass
# def air_enthalpy(T,conc_mass,Mm_a): #==> a chequer si c est pas diviser par Mm_a ou divisé par molar_mass
#     enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
#     molar_mass = np.array([0.028,0.044,0.018,0.032])
#     h_air = sum(conc_mass*enthalpies);#kJ/mol
#     return h_air/Mm_a #kJ/kg
# def air_enthalpy(T,conc_mass,Mm_a): #==> a chequer si c est pas diviser par Mm_a ou divisé par molar_mass
#     enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
#     h_air = sum(conc_mass*enthalpies);#kJ/mol
#     return h_air/Mm_a #kJ/kg
# def air_entropy(T,conc_mass,Mm_a):
#     entropies = np.array([N2.S(T),CO2.S(T),H2O.S(T),O2.S(T)])
#     S_air = sum(conc_mass*entropies);#J/mol/K
#     return S_air/Mm_a #kJ/kg

def cp_air_T(T,conc_mass,Mm_a):#J/kg/K
    return cp_air(T,conc_mass,Mm_a)/T;

def janaf_integrate_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)*dt)
def cp_mean_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)/len(values)) #  cp_mean [J/kg/K]
def janaf_integrate(f,T1,T2,dt): #==> pour calculer enthalpie
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt) # int(cp)dt [J/mol/K]]


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
    return I,cp_mean,heat_const
