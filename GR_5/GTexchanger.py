from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GT_arguments as GT_arg;
import GTcomb_arguments as GTcomb_arg
import combustionGT as comb;
import exchanger
import useful

"""
Convention :  All temperatures in input are given in °C but are immediatly changed into K
inside of the function. output temperatures are also in °C
"""

def GT(GT_input):
    """
     GT Gas turbine modelisation
     GT(P_e,options,display) compute the thermodynamics states for a Gas
     turbine based on several inputs (given in OPTION) and based on a given
     electricity production P_e. It returns the main results. It can as well
     plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)

     In this GT, there is a heat exchanger between the compressor and the
     combustion chamber. The incoming air is heated up using the hot gasses
     of the exhaust. In this model we choosed to heat up the air after the
     compressor 100 K (T2r = T2+100)

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
     OUTPUTS : outputs are specified in the class GT_outputs in GT_arguments.py
    """
    arg_in = GT_input;

    ## Check input arguments
    # ======================
    Pe = arg_in.Pe;
    if Pe ==-1.:
        Pe = 50e3;#50MWe
    T_ext = arg_in.T_ext;
    if T_ext ==-1.:
        T_ext = 15;#15°C
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
    if T0 ==-1.:
        T0 = 15 #[°C]
    k_mec = arg_in.k_mec
    if k_mec==-1.:
        k_mec = 0
    kcc = arg_in.k_cc
    if kcc == -1.:
        kcc = 1

    T3 =T3 +273.15 #[K]
    T_ext = T_ext +273.15 #[K]
    T0 = T0+273.15 #[K]
    """
    0) Air data
    """
    Mm_O2 = 0.032;#kg/mol
    Mm_N2 = 0.028;#kg/mol
    conc_O2 = 0.21;# 21% in molar
    conc_N2 = 0.79;# 79% in molar
    Cp_a,gamma,R= useful.air_mixture(T0)
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    conc_mass1=np.array([conc_N2*Mm_N2/Mm_a,0,0,conc_O2*Mm_O2/Mm_a])
    Ra = 8.31/Mm_a
    #polytropic coefficient (compressor :m>gamma , turbine : gamma >m) first estimation
    m_t = (-eta_pit*(gamma-1)/gamma+1)**(-1)
    m_c = (1-1/eta_pic*(gamma-1)/gamma)**(-1)

    """
    1) compressor

    """
    dt = 0.01
    T1=T_ext
    p1 = 1.0 #bar
    h1 = useful.janaf_integrate_air(useful.cp_air,conc_mass1,Mm_a,T0,T1,dt)/1000 #kJ/kg
    s1 = useful.janaf_integrate_air(useful.cp_air_T,conc_mass1,Mm_a,T0,T1,dt)#J/kg/K
    e1 = h1-T0*s1/1000 #kJ/kg_in

    p2 = r*p1

    # calculation of T2 via iterations
    T2 = T1*(r)**((m_c-1)/m_c) #first estimation
    iter = 1
    error = 1
    dt = 0.1

    while iter < 50 and error >0.01 :#pg119
        exposant_c=1/eta_pic*(Ra)/useful.cp_mean_air(useful.cp_air,conc_mass1,Mm_a,T1,T2,0.001)  # na-1/na :  formule du livre (3.22)
        T2_new = T1*r**(exposant_c)
        iter=iter+1
        error = abs(T2_new-T2)
        T2=T2_new

    s2 = useful.janaf_integrate_air(useful.cp_air_T,conc_mass1,Mm_a,T0,T2,0.001)-Ra*np.log(r)
    h2 = useful.janaf_integrate_air(useful.cp_air,conc_mass1,Mm_a,T0,T2,dt)/1000
    e2 = h2-T0*s2/1000

    deltah_c = h2-h1 #kJ/kg
    deltas_c1 = s2-s1
    delta_ex_c = e2-e1

    """
    1.b) Preheating
    """
    T2r = T2 +100 # we choose to heat it up 100 K
    h2r = useful.janaf_integrate_air(useful.cp_air,conc_mass1,Mm_a,T0,T2r,dt)/1000
    s2r = useful.janaf_integrate_air(useful.cp_air_T,conc_mass1,Mm_a,T0,T2r,0.001)-Ra*np.log(r)#J/kg/K
    e2r =  h2-T0*s2/1000#kJ/kg_air

    """
     2 ) Combustion
    """
    p3 = p2*kcc

    comb_inputs = GTcomb_arg.comb_input(h_in=h2r,T_in = T2r-273.15,inversion=True,T_out=T3-273.15, k_cc = kcc, r=r )
    comb_outputs = comb.combustionGT(comb_inputs)

    T3=comb_outputs.T_out+273.15 #[K]
    lambda_comb = comb_outputs.lambda_comb
    ma1 = comb_outputs.ma1

    #new molar mass and concentration of the exhaust gasses
    Mm_af = comb_outputs.Mm_af
    Rf = comb_outputs.R_f
    conc_mass2 = np.array([comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of,comb_outputs.m_O2f])
    massflow_coefficient = 1+1/(ma1*lambda_comb) #kg_fu/kg_air

    h3 = useful.janaf_integrate_air(useful.cp_air,conc_mass2,Mm_af,T0,T3,0.001)/1000#kJ/kg
    s3 = useful.janaf_integrate_air(useful.cp_air_T,conc_mass2,Mm_af,T0,T3,0.001)-Rf*np.log(kcc*r) #J/K/kg
    e3 = h3-T0*s3/1000 #kJ/kg_in
    delta_exer_comb = massflow_coefficient*e3-e2r #kJ/kg_air

    """
    3)  Turbine
    """
    p4 = p3/(r*kcc)

    # calculation of T4 via iterations
    T4 = T3*(1/(r*kcc))**((m_t-1)/m_t)
    iter = 1
    error = 1
    dt = 0.1

    while iter < 50 and error >0.01 :
        exposant_t=eta_pit*(Rf)/useful.cp_mean_air(useful.cp_air,conc_mass2,Mm_af,T4,T3,0.01)  # na-1/na :  formule du livre (3.25)
        T4_new = T3*(1/(kcc*r))**(exposant_t)
        iter=iter+1
        error = abs(T4_new-T4)
        T4=T4_new

    h4 = useful.janaf_integrate_air(useful.cp_air,conc_mass2,Mm_af,T0,T4,0.001)/1000 #kJ/kg_f
    deltah_t = h4-h3 #<0# kJ/kg_f
    s4 = useful.janaf_integrate_air(useful.cp_air_T,conc_mass2,Mm_af,T0,T4,0.001) #J/K/kg
    e4 = h4-T0*s4/1000# kJ/kg_f
    delta_exer_t = e4-e3# kJ/kg_f
    deltas_t = s4-s3# kJ/kg_f

    """
    4) Mechanical work (travail moteur)
    """
    # Mechanical work
    Wm = -(deltah_c+(1+1/(lambda_comb*ma1))*deltah_t) #kJ/kg_in
    #heat in
    Q_comb = massflow_coefficient*h3-h2r #kJ/kg_in

    """
    5) Energetic efficiencies and mass flux
    """
    # efficiencies
    eta_cyclen  =Wm/Q_comb
    eta_mec =1-k_mec* (massflow_coefficient*abs(deltah_t)+deltah_c)/(massflow_coefficient*abs(deltah_t)-deltah_c) # Pe/Pm = 1-k_mec*(Pmt+Pmc)/(Pmt-Pmc)
    eta_toten = eta_cyclen*eta_mec

    #massflow
    mf_in = Pe/(Wm*eta_mec)#kg/s
    mf_out = mf_in*massflow_coefficient #kg/s
    mf_c = mf_out-mf_in #kg/s

    """
    3.b) Heat exhanger
    """
    exchanger_inputs = GTcomb_arg.exchanger_input(Mflow_air_in = mf_in,Mflow_f_in = mf_out,T_air_in=T2-273.15,
                                                        T_f_in = T4-273.15, comb_lambda=lambda_comb,U = 3)
    exchanger_outputs = exchanger.heatexchanger(exchanger_inputs,T2r-273.15)

    T5 = exchanger_outputs.T_f_out+273.15
    h5 = useful.janaf_integrate_air(useful.cp_air,conc_mass2,Mm_af,T0,T5,0.001)/1000 #kJ/kg
    s5 = useful.janaf_integrate_air(useful.cp_air_T,conc_mass2,Mm_af,T0,T5,0.001) #J/K/kg
    e5 = h5-T0*s5/1000
    Q_pre = exchanger_outputs.Q #kJ

    """
    5) Calculation of the energetic fluxes
    """
    P_in = h1*mf_in #[kW]
    P_c = deltah_c*mf_in #[kW]
    P_pre = Q_pre*mf_in #[kW]
    P_comb = Q_comb*mf_in #[kW]
    P_t = -deltah_t*mf_in*massflow_coefficient
    P_out = h5*massflow_coefficient*mf_in #[kW]
    P_fmec = P_t-P_c-Pe
    Pm = P_t-P_c

    """
    7) Calculation of the losses (compressor, comb, turbine, exhaust)
    """
    P_ech = P_comb-Pm
    #compressor losses
    L_c = mf_in*T0*deltas_c1/1000 #[kW]
    # combustion losses
    ec = comb_outputs.e_c
    L_comb = mf_in*e2r-mf_out*e3+mf_c*ec #[kW]
    #turbine losses
    L_t = mf_out*T0*deltas_t/1000 #[kW]
    #exhaust losses
    L_exhaust = mf_out*e5 #-mf_in*e1
    #exchange losses :
    L_exchanger = e2*mf_in+e4*mf_out-e5*mf_out-e2r*mf_in

    """
    8) Exergetic efficiencies
    """
    eta_cyclex = Pm/(mf_out*e3-mf_in*e2r)
    eta_totex = Pe/(mf_out*h3-mf_out*h2r)
    eta_rotex = Pm/(mf_out*(e3-e4)-mf_in*(e2-e1))
    eta_combex = (mf_out*e3-mf_in*e2r)/(mf_c*ec)
    eta_cex = delta_ex_c/deltah_c
    eta_dex = deltah_t/delta_exer_t
    """
    9) Define output arguments
    """
    outputs = GT_arg.GT_outputs();
    outputs.eta[0] = eta_cyclen;
    outputs.eta[1] = eta_toten;
    outputs.eta[2] = eta_cyclex;
    outputs.eta[3] = eta_totex;
    outputs.eta[4] = eta_rotex;
    outputs.eta[5] = eta_combex;
    outputs.daten[0] = P_ech; #[kW]
    outputs.daten[1] = P_fmec;#[kW]
    outputs.dat[0:]= [[T1-273.15,T2-273.15,T3-273.15,T4-273.15],[p1,p2,p3,p4],[h1,h2,h3,h4],[s1/1000,s2/1000,s3/1000,s4/1000],[e1,e2,e3,e4]]
    outputs.massflow[0:] = [mf_in,mf_c,mf_out]
    outputs.combustion.fum[0:]=np.array([comb_outputs.m_O2f,comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of])*mf_out
    outputs.combustion.Lambda = lambda_comb
    outputs.combustion.LHV = comb_inputs.LHV
    outputs.combustion.e_c = comb_outputs.e_c
    outputs.combustion.Cp_g = useful.cp_air(400,conc_mass2,Mm_af)/1000

    """
    10) Pie charts and cycle graphs
    """

    # pie chart of the energie flux in the cycle
    fig,ax =  plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    data = [Pe,P_fmec,P_out]
    labels = ['Useful power {v} [MW]'.format(v=round(Pe/1000)),'Mechanical losses {v} [MW]'.format(v=round(P_fmec/1000)),'Exhaust losses {v} [MW]'.format(v=round(P_out/1000))]

    ax.pie(data,labels = labels,autopct='%1.2f%%',startangle = 90)
    ax.set_title("Primary energetic flux "+ str(round(P_comb/10**3)) + "[MW]")

    # pie chart of the exergie flux in the cycle
    fig2,ax =  plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    data = [Pe,P_fmec,L_t+L_c,L_exhaust,L_exchanger,L_comb]
    labels = ['Useful power {v} [MW]'.format(v=round(Pe/1000)),'Mechanical losses {v} [MW]'.format(v=round(P_fmec/1000)),'\n Compressor and turbine losses {v} [MW]'.format(v=round((L_t+L_c)/1000)),
             'Exhaust losses {v} [MW]'.format(v=round(L_exhaust/1000)),'Heat exchanger losses {v} [MW]'.format(v=round(L_exchanger/1000)), 'Combustion losses {v} [MW]'.format(v=round(L_comb/1000))]


    ax.pie(data,labels = labels,autopct="%1.2f%%",startangle = 90)
    ax.set_title("Primary exergetic flux "+ str(round(ec*mf_c/10**3)) + "[MW]")

    # T S graph of the cycle
    Ta = np.linspace(T1,T2,50)
    Tc = np.linspace(T4,T3,50)
    Sa= np.zeros(len(Ta))
    Sc = np.zeros(len(Tc))

    for i in range(len(Ta)):
        Sa[i] = s1+(1-eta_pic)*useful.janaf_integrate_air(useful.cp_air_T,conc_mass1,Mm_a,T1,Ta[i],dt)
        Sc[i] = s4-(1-eta_pit)/eta_pit*useful.janaf_integrate_air(useful.cp_air_T,conc_mass2,Mm_af,T4,Tc[i],dt)

    Sb=np.linspace(s2r,s3,50)
    Sb_pre = np.linspace(Sa[-1],s2r,40)
    a,b= np.polyfit([s2r,s3],[T2r,T3],1)
    a_pre,b_pre= np.polyfit([Sa[-1],s2r],[T2,T2r],1)
    Sd=np.linspace(s1,s5,50)
    a2,b2= np.polyfit([s1,s5],[T1,T5],1)
    Sd_pre=np.linspace(s5,s4,40)
    a2_pre,b2_pre= np.polyfit([s5,s4],[T5,T4],1)

    fig3,ax1 = plt.subplots()
    ax1.plot(Sa,Ta-273.15,Sc,Tc-273.15,Sb,a*Sb+b-273.15,Sd,a2*Sd+b2-273.15,Sb_pre,a_pre*Sb_pre+b_pre-273.15,Sd_pre,a2_pre*Sd_pre+b2_pre-273.15)
    ax1.scatter([s1,Sa[-1],s2r,s3,s4,s5],[T1-273.15,T2-273.15,T2r-273.15,T3-273.15,T4-273.15,T5-273.15],s=10,label='States')
    ax1.set_xlabel('Entropy [J/kg/K]')
    ax1.set_ylabel('Tempearature [°C]')
    ax1.grid(True)
    ax1.legend()
    ax1.set_title('T S graph of the gaz turbine cycle')

    # p v graph of the cycle
    pa = np.linspace(p1,p2,50)
    pb = np.linspace(p4,p3,50)
    va = R/pa*(pa/p1)**(exposant_c)*T1
    vb = R/pb*(pb/p4)**(exposant_t)*T4
    m,b = np.polyfit([va[-1],vb[-1]], [pa[-1],pb[-1]], 1)
    vc = np.linspace(va[-1],vb[-1],10)
    pc = m*vc+b
    vd = np.linspace(va[0],vb[0],10)

    fig4,ax2=plt.subplots()
    ax2.plot(va,pa,vb,pb,vc,pc,vd,p4*np.ones(len(vd)))
    ax2.scatter([R*T1/p1,R*T2/p2,R*T2r/p2,R*T3/p3,R*T4/p4,R*T5/p4],[p1,p2,p2,p3,p4,p4],s=10,label='States')
    ax2.set_xlabel('Specific volume $[m^3/kg]$')
    ax2.set_ylabel('Pressure [bar]')
    ax2.grid(True)
    ax2.legend()

    fig = [fig,fig2,fig3,fig4]
    outputs.fig = fig
    if (GT_input.Display == 1):
        plt.show()

    return outputs
# GT_simple_outputs = GT(GT_arg.GT_input(Pe = 50e3,k_mec =0.015, T_ext=15,T_0 = 15,r=18.,k_cc=0.95,T3 = 1400,Display =1));
# print(GT_simple_outputs.dat)
