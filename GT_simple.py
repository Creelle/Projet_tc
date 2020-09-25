
from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;

import GT_arguments as GT_arg;
import combustionGT as comb;

O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');

def air_mixture(T):#kJ/kg/K
    Mm_O2 = 0.032;#kg/mol
    Mm_N2 = 0.028;#kg/mol
    conc_O2 = 0.21;# 21% in molar
    conc_N2 = 0.79;# 79% in molar
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg/K
    R = 8.31/Mm_a/1000
    gamma = Cp/(Cp-R)
    return Cp,gamma ;

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
    T2 = T_ext*(r)**((m_c-1)/m_c)
    deltah_c = Cp_a*(T2-T_ext)
    print(T2,deltah_c)
    comb.combustionGT(GT_arg.comb_input())

    ## define output arguments
    # ======================
    outputs = GT_arg.GT_outputs();
    outputs.eta[1] = 0.35;

    # Your job

    return outputs;




#tests
GT_simple_outputs = GT_simple(GT_arg.GT_input());
