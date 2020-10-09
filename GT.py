
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

def air_mixture(T):#kJ/kg/K
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg/K
    R = 8.31/Mm_a/1000
    gamma = Cp/(Cp-R)
    return Cp,gamma;


#fonction qui donne l enthalpie (kJ/kg), T temperature concentration massique mass : array (N2 , CO2, H20,O2)
#molar mass
def air_enthalpy(T,conc_mass,Mm_a):
    enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
    h_air = sum(conc_mass*enthalpies);#kJ/mol/K
    return h_air/Mm_a #kJ/kg

def air_entropy(T):
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    entropy_air = m_O2 * O2.S(T) + N2.S(T) * m_N2;#kJ/mol/K
    return entropy_air/Mm_a #kJ/kg/K

def janaf_integrate(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt)

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
        a faire:
        merge combustionGT avec simpleGT
        exergie pour simpleGT
        optimiser le taux de compression jouant sur le taux de compression et le lambda
        preheating ==> modelisation d un echangeur
        (humidity chequ)
        pychart
        des graphes T s et pv des etats dans la turbine
        (faire plusieurs etages de compression et analyse au niveau exergetique et energetique pour avoir
        si ca change quelque chose)
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
    k_mec = arg_in.k_mec


    ## preliminary data (air) ==> find gamma
    # ======================
    # cp air at 15°C (298K): [kJ/mol/K]
    Cp_a,gamma= air_mixture(T0+273.15)
    #molar mass entry
    Mm_a = Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    conc_mass1=np.array([conc_N2,0,0,conc_O2])
    #coeff polytroique (compresseur :m>gamma , turbine : gamma >m)
    gamma = Cp_a/(Cp_a-287.1)
    m_t = (-eta_pit*(gamma-1)/gamma+1)**(-1)
    m_c = (1-1/eta_pic*(gamma-1)/gamma)**(-1)
    # cycle definition
    # =================
    #1) compressor
    #exergie = dh-T0 * ds
    T1=T_ext
    p1 = 1.0 #bar
    h1 = air_enthalpy(T1,conc_mass1,Mm_a)
    s1 = air_entropy(T1)

    p2 = r*p1
    T2 = T_ext*(r)**((m_c-1)/m_c)
    s2 = air_entropy(T2)

    h2 = air_enthalpy(T2,conc_mass1,Mm_a)
    deltah_c = h2-h1
    #deltah_c = Cp_a*(T2-T1) # delta_h =  w_m compression
    #deltas_c = Cp_a*np.log(T2/T1)*1000

    #delats_c = janaf_integrate(air_mixture2,T1,T2,0.0001) ca ne marche pas


    # 2) combustion
    p3 = p2*k_cc
    comb_outputs = comb.combustionGT(GT_arg.comb_input(h_in=h2,T_in = T2,lambda_comb = 5))
    T3=comb_outputs.T_out
    lambda_comb = comb_outputs.lambda_comb
    ma1 = comb_outputs.ma1
    Mm_af = comb_outputs.Mm_af
    conc_mass2 = np.array([comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of,comb_outputs.m_O2f])
    h3 = air_enthalpy(T3,conc_mass2,Mm_af) #kJ/kg
    massflow_coefficient = 1+1/(ma1*lambda_comb)
    Q=h3-h2
    s3 = air_entropy(T3)

    # 3)  combustion
    p4 = p3/r
    T4 = T3*(1/(r*kcc))**((m_t-1)/m_t)

    h4 = air_enthalpy(T4,conc_mass2,Mm_af)
    deltah_t = h4-h3
    print(Cp_a*(h4-h3),deltah_t)
    s4 = air_entropy(T4)

    #travail moteur
    Wm = -(deltah_c+(1+1/(lambda_comb*ma1))*deltah_t) #kJ/kg
    #apport calorifique
    Q_comb = (1+1/(lambda_comb*ma1))*(h3-h2)
    # autre variable utile : X= (p2/p1)**((gamma-1)/gamma))
    print('1',deltah_c+deltah_t,h4-h3,h2-h1, s3-s4)

    ##====================
    # calcul des rendements
    eta_cyclen  =Wm/Q_comb
    eta_mec = 1-k_mec#(Wm - k_mec)/Wm
    eta_gen = 1
    eta_toten = eta_cyclen*eta_mec*eta_gen #
    print('2',eta_cyclen, (deltah_t+deltah_c)/Q)
    #massflow calcul # on neglige m  flow combustion
    mf_in = Pe/(Wm*eta_mec)#
    mf_out = mf_in*massflow_coefficient
    ## define output arguments
    # ======================
    outputs = GT_arg.GT_outputs();
    outputs.eta[0] = eta_cyclen;
    outputs.eta[1] = eta_toten;
    outputs.dat[0:4]= [[T_ext,T2,T3,T4],[p1,p2,p3,p4],[h1,h2,h3,h4],[s1,s2,s3,s4]]
    outputs.massflow[0:] = [mf_in,mf_out-mf_in,mf_out]
    outputs.combustion.fum[0:]=np.array([comb_outputs.m_O2f,comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of])*mf_out
    # Your job
    print('3',outputs.massflow)
    return outputs;




#tests
GT_simple_outputs = GT_simple(GT_arg.GT_input());
