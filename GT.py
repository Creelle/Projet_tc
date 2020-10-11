
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

def cp_air(T,conc_mass,Mm_a):
    cps = np.array([N2.cp(T),CO2.cp(T),H2O.cp(T),O2.cp(T)])
    molar_mass = np.array([0.028,0.044,0.018,0.032])
    cp_air = np.dot(conc_mass/molar_mass,cps);#J/mol/K
    return cp_air#/Mm_a #J/kg


#fonction qui donne l enthalpie (kJ/kg), T temperature concentration massique mass : array (N2 , CO2, H20,O2)
#molar mass
def air_enthalpy(T,conc_mass,Mm_a): #==> a chequer si c est pas diviser par Mm_a ou divisé par molar_mass
    enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
    molar_mass = np.array([0.028,0.044,0.018,0.032])
    h_air = sum(conc_mass/molar_mass*enthalpies);#kJ/mol
    return h_air#/Mm_a #kJ/kg

def air_entropy(T,conc_mass,Mm_a):
    entropies = np.array([N2.S(T),CO2.S(T),H2O.S(T),O2.S(T)])
    S_air = sum(conc_mass*entropies);#J/mol/K
    return S_air/Mm_a #kJ/kg
def cp_air_T(T,conc_mass,Mm_a):#J/kg/K

    return cp_air(T,conc_mass,Mm_a)/T;

def exergie_air(T,conc_mass,Mm_a):
    T0=288.15
    #molar_mass = np.array([0.028,0.044,0.018,0.032])
    enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
    entropies = np.array([N2.S(T),CO2.S(T),H2O.S(T),O2.S(T)])
    enthalpies0 = np.array([N2.hef(T0),CO2.hef(T0),H2O.hef(T0),O2.hef(T0)])
    entropies0 = np.array([N2.S(T0),CO2.S(T0),H2O.S(T0),O2.S(T0)])
    exergies = (enthalpies-enthalpies0)*1000-T0*(entropies-entropies0) #J/mol
    e_air = sum(conc_mass*exergies)/1000/Mm_a #kJ/kg
    return e_air

def janaf_integrate_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)*dt)
def cp_mean_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)/len(values)) #  cp_mean [J/kg/K]
def janaf_integrate(f,T1,T2,dt): #==> pour calculer enthalpie
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt) # int(cp)dt [J/mol/K]]


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

        a faire dans l immediat
        changer la formule de l entropy_air
        janaf integrate air
        formule d exergetique
        rendements exergetique
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
    kcc = arg_in.k_cc

    """
    ## preliminary data (air) ==> find gamma
    # ======================
    # cp air at 15°C (298K): [kJ/mol/K]
    """
    Cp_a,gamma= air_mixture(T0)
    Mm_a = Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    conc_mass1=np.array([conc_N2*Mm_N2/Mm_a,0,0,conc_O2*Mm_O2/Mm_a])
    Ra = 8.31/Mm_a
    #coeff polytroique (compresseur :m>gamma , turbine : gamma >m) en premiere estimation
    m_t = (-eta_pit*(gamma-1)/gamma+1)**(-1)
    m_c = (1-1/eta_pic*(gamma-1)/gamma)**(-1)
    #on va recalculer m_c et m_t en utilisant la definition du livre page 118 (3.20) (3.25 en faisant des iterations)


    """
    1) compressor

    """
    T1=T_ext # a changer lors du preaheating
    p1 = 1.0 #bar
    h1 = air_enthalpy(T1,conc_mass1,Mm_a)
    s1 = air_entropy(T1,conc_mass1,Mm_a)

    p2 = r*p1

    # calcul de T2 par iteration
    T2 = T1*(r)**((m_c-1)/m_c) #premiere estimation
    iter = 1
    error = 1
    dt = 0.1

    while iter < 50 and error >0.01 :#pg119
        exposant_c=1/eta_pic*(Ra)/cp_mean_air(cp_air,conc_mass1,Mm_a,T1,T2,0.01)  # na-1/na :  formule du livre (3.2)
        T2_new = T1*r**(exposant_c)
        iter=iter+1
        error = abs(T2_new-T2)
        T2=T2_new

    s2 = air_entropy(T2,conc_mass1,Mm_a)
    h2 = air_enthalpy(T2,conc_mass1,Mm_a)
    deltah_c = h2-h1 #kJ/kg
    deltah_c2 = janaf_integrate_air(cp_air,conc_mass1,Mm_a,T1,T2,0.01)
    print('enthalpy comparaison',deltah_c,deltah_c2)

    deltas_c1 = s2-s1
    deltas_c2 = Cp_a*np.log(T2/T1)*1000 #J/K/kg
    deltas_c3 = janaf_integrate_air(cp_air_T,conc_mass1,Mm_a,T1,T2,0.01)
    print('entropy comparaison',deltas_c1,deltas_c2,deltas_c3)

    #deltah_c = Cp_a*(T2-T1) # delta_h =  w_m compression
    #deltas_c = Cp_a*np.log(T2/T1)*1000


    """
     2 ) combustion
    """

    p3 = p2*kcc

    comb_outputs = comb.combustionGT(GT_arg.comb_input(h_in=h2,T_in = T2,lambda_comb =1.65 ))
    T3=comb_outputs.T_out
    lambda_comb = comb_outputs.lambda_comb
    ma1 = comb_outputs.ma1
    Mm_af = comb_outputs.Mm_af
    Rf = comb_outputs.R_f
    conc_mass2 = np.array([comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of,comb_outputs.m_O2f])

    h3 = air_enthalpy(T3,conc_mass2,Mm_af) #kJ/kg
    massflow_coefficient = 1+1/(ma1*lambda_comb)
    s3 = air_entropy(T3,conc_mass2,Mm_af)
    print('s3',s3, s2+janaf_integrate_air(cp_air_T,conc_mass2,Mm_af,T2,T3,0.01))
    """
    3)  detente
    """
    p4 = p3/(r*kcc)


    # calcul de T4 par iteration
    T4 = T3*(1/(r*kcc))**((m_t-1)/m_t)
    print(T3)
    print(T4)
    iter = 1
    error = 1
    dt = 0.1

    while iter < 50 and error >0.01 :#pg119
        exposant_t=eta_pit*(Rf)/cp_mean_air(cp_air,conc_mass2,Mm_af,T4,T3,0.01)  # na-1/na :  formule du livre (3.2)
        T4_new = T3*(1/(kcc*r))**(exposant_t)
        iter=iter+1
        error = abs(T4_new-T4)
        T4=T4_new
        print(T4)


    h4 = air_enthalpy(T4,conc_mass2,Mm_af)
    deltah_t = h4-h3
    s4 = air_entropy(T4,conc_mass2,Mm_af)

    """
    4) travail moteur et rendements
    """
    #travail moteur
    Wm = -(deltah_c+(1+1/(lambda_comb*ma1))*deltah_t) #kJ/kg_in
    #apport calorifique
    Q_comb = massflow_coefficient*h3-h2
    # autre variable utile : X= (p2/p1)**((gamma-1)/gamma))

    ##====================
    # calcul des rendements
    eta_cyclen  =Wm/Q_comb
    eta_mec = 1-k_mec#(Wm - k_mec)/Wm
    eta_gen = 1
    eta_toten = eta_cyclen*eta_mec*eta_gen #
    print('2',eta_cyclen, (deltah_t+deltah_c)/Q_comb)
    #massflow calcul # on neglige m  flow combustion
    mf_in = Pe/(Wm*eta_mec)#
    mf_out = mf_in*massflow_coefficient

    """
    5) define the exergy
    """
    e1 = exergie_air(T1,conc_mass1,Mm_a)
    e2 = exergie_air(T2,conc_mass1,Mm_a)
    e3 = exergie_air(T3,conc_mass2,Mm_af)
    e4 = exergie_air(T4,conc_mass2,Mm_af)
    """
    last) define output arguments
    """
    outputs = GT_arg.GT_outputs();
    outputs.eta[0] = eta_cyclen;
    outputs.eta[1] = eta_toten;
    outputs.dat[0:]= [[T1,T2,T3,T4],[p1,p2,p3,p4],[h1,h2,h3,h4],[s1,s2,s3,s4],[e1,e2,e3,e4]]
    outputs.massflow[0:] = [mf_in,mf_out-mf_in,mf_out]
    outputs.combustion.fum[0:]=np.array([comb_outputs.m_O2f,comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of])*mf_out

    return outputs;
"""
attention, la temperature de reference dans janaf n est pas 288.15 mais 298.15
"""

GT_simple_outputs = GT_simple(GT_arg.GT_input(Pe = 230e3,T_ext=288.15,r=18.));
print(GT_simple_outputs.dat)
print(GT_simple_outputs.massflow)
