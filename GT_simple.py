
from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;

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

def air_mixture(T):#kJ/kg/K
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg/K
    R = 8.31/Mm_a/1000
    gamma = Cp/(Cp-R)
    return Cp,gamma ;

#fonction qui donne l enthalpie (kJ/kg)
def air_enthalpy(T):
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    h_air = m_O2 * O2.hef(T) + N2.hef(T) * m_N2;#J/mol/K
    return h_air/Mm_a #kJ/kg

def air_entropy(T):
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    entropy_air = m_O2 * O2.S(T) + N2.S(T) * m_N2;#J/mol/K
    return entropy_air/Mm_a #kJ/kg/K

def GT_simple(GT_input):
    """
     GT Gas turbine modelisation
     GT(P_e,options,display) compute the thermodynamics states for a Gas
     turbine based on several inputs (given in OPTION) and based on a given
     electricity production P_e. It returns the main results. It can as well
     plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)

     INPUTS (some inputs can be dependent on others => only one of these 2 can
             be activated) Refer to Fig 3.1 from reference book (in english)
     P_E = electrical power output target [kW]
     OPTIONS is a structure containing :
       -options.T_ext [°C] : External temperature
       -options.r     [-] : Comperssion ratio
                            chamber
       -options.T_3   [°C] : Temperature after combustion (before turbine)
       -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for compression
       -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for expansion
    """
    arg_in = GT_input;

    ## Check input arguments
    # ======================
    Pe = arg_in.Pe;
    if Pe ==-1.:
        Pe = 50e3;#50MWe
    T_ext = arg_in.T_ext;
    if T_ext ==-1.:
        T_ext = 288.15;#15°C
    r = arg_in.r;
    if r ==-1.:
        r = 10;#compression ratio = 10;
    T3 =arg_in.T3;
    eta_pic = arg_in.eta_PiC;
    if eta_pic ==-1.:
        eta_pic = 0.9;#max temperature = 1050°C
    eta_pit = arg_in.eta_PiT;
    if eta_pit ==-1.:
        eta_pit = 0.9;#max temperature = 1050°C
    T0=arg_in.T_0


    ## preliminary data (air) ==> find gamma
    # ======================
    # cp air at 15°C (298K): [kJ/mol/K]
    Cp_a,gamma = air_mixture(T0+273.15)
    #coeff polytroique (compresseur :m>gamma , turbine : gamma >m)
    m_t = (-eta_pit*(gamma-1)/gamma+1)**(-1)
    m_c = (1-1/eta_pic*(gamma-1)/gamma)**(-1)
    ## cycle definition
    # =================
    #1) compressor
    T1=T_ext
    p1 = 1.0 #bar
    h1 = air_enthalpy(T1)
    s1 = air_entropy(T1)

    p2 = r*p1
    T2 = T_ext*(r)**((m_c-1)/m_c)
    s2 = air_entropy(T2)
    #h2 = air_enthalpy(T2)
    deltah_c = Cp_a*(T2-T1) # delta_h =  w_m compression
    deltas_c = Cp_a*np.log(T2/T1)*1000
    h2 = h1+deltah_c
    # 2) combustion
    p3 = p2
    h3 = air_enthalpy(T3)
    Q=h3-h2
    s3 = air_entropy(T3)

    # 3)  combustion
    p4 = p3/r
    T4 = T3*(1/r)**((m_t-1)/m_t)
    deltah_t = Cp_a*(T4-T3)
    #h4 = air_enthalpy(T4)
    h4 = h3+deltah_t
    s4 = air_entropy(T4)
    # autre variable utile : X= (p2/p1)**()(gamma-1)/gamma)
    print(deltah_c+deltah_t,h4-h3,h2-h1, s3-s4)
    comb.combustionGT(GT_arg.comb_input())
    ##====================
    # calcul des rendements
    eta_cyclen  = 1-(h4-h1)/(h3-h2)
    print(eta_cyclen, (deltah_t+deltah_c)/Q)
    ## define output arguments
    # ======================
    outputs = GT_arg.GT_outputs();
    outputs.eta[1] = 0.35;

    # Your job

    return outputs;




#tests
GT_simple_outputs = GT_simple(GT_arg.GT_input());
