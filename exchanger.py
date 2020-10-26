from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GTcomb_arguments as GT_comb_arg;
import combustionGT as comb;

db = janaf.Janafdb();
Mm_CH4 = 0.016; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044
O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
CH4=db.getphasedata('CH4',phase ='g');
CO = db.getphasedata('CO',phase ='g')

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



def heatexchanger(exchanger_input,T_air_out):
    """
    On connait le mass flow d'air in et le mass flow de fumée in ainsi que tous les cp correspondants
    --> on les connait car on veut une certaine quantité of power à la sortie directement lié au massflow dans la GT.
    On connait la température in de l'air et la température in des fumées (calculée par combustionGT).
    J'indique dans la défintion de la fonction, la température out de l'air que je souhaite à la sortie de l'échangeur
    Return la température des fumées à la sortie de l'échangeur ainsi que la surface d'échange nécessaire (donné grâce à Hausbrand)

    ATTENTION LE CODE A UNE LIMITE : T_f_in DOIT ÊTRE INFERIEUR À 988K si je veux T_air_out à 500K.
    Si j'augmente T_air_out alors je peux augmenter T_f_in
    Est-ce un problème numérique ou physique?
    """
    U = exchanger_input.U
    Mflow_air_in = exchanger_input.Mflow_air_in
    Mflow_f_in = exchanger_input.Mflow_f_in
    T_air_in = exchanger_input.T_air_in +273.15 #[K]
    T_f_in = exchanger_input.T_f_in +273.15 #[K]
    lambda_comb = exchanger_input.comb_lambda

    x_O2a = 0.21
    x_N2a = 0.79
    coeff = x_N2a/x_O2a
    T0 =  288.15 #[K]

    Mm_a = x_O2a * Mm_O2 + x_N2a * Mm_N2 # [kg/mol_air]
    ma1 =  Mm_a/Mm_CH4 * 2/x_O2a # kg_air/kg_CH4 = proportion d air entrant vs combustible

    # A la sortie :
    coeff_stochio = np.array([2*lambda_comb*coeff,1,2,2*(lambda_comb-1)]) # N2- CO2 - H2O - O2
    total_n = sum(coeff_stochio) # nombre de moles total
    molar_conc = coeff_stochio/total_n # concentration des elements mol_co2/mol_t
    molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
    Mm_af = sum(molar_conc*molar_mass) #somme ponderé des masse molaire
    mass_conc = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot


    dt = 0.01
    iter = 1
    error = 20

    T_air_out = T_air_out +273.15 # K
    Q = Mflow_air_in * janaf_integrate(cp_air,T_air_in,T_air_out,dt)
    print("Q : ",Q,'[J]')

    #je fais une première estimation de T_f_out pour la formule itérative ci-dessous
    #cette première estimation se fait à cp constant pour l'intégration
    cps_out = np.array([cpN2(288.15),cpCO2(288.15),cpH2O(288.15),cpO2(288.15)])#premiere estimation
    cp_f = np.dot(cps_out,mass_conc)/Mm_af #J/kg_fumée
    T_f_out = (-Q/(Mflow_f_in*cp_f)) + T_f_in
    T_f_out_secours = T_f_out

    #Nous avons itéré en utilisant l'interpolation par un polynome du troisième degré des différents cp des composants
    #Les autres méthodes d'intégration ne nous donnaient pas une précision assez élevée.

    while iter<50 and error>0.01 :
        iter = iter+1
        heat_const_N2 = cp_Iconstants('N2',T_f_out,T_f_in)
        heat_const_CO2 = cp_Iconstants('CO2',T_f_out,T_f_in)
        heat_const_H2O = cp_Iconstants('H2O',T_f_out,T_f_in)
        heat_const_O2 = cp_Iconstants('O2',T_f_out,T_f_in)
        int_N2 = (1/2)*heat_const_N2[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_N2[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_N2[2][0]*(T_f_in**4-T_f_out**4)
        int_CO2 = (1/2)*heat_const_CO2[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_CO2[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_CO2[2][0]*(T_f_in**4-T_f_out**4)
        int_H2O = (1/2)*heat_const_H2O[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_H2O[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_H2O[2][0]*(T_f_in**4-T_f_out**4)
        int_O2 = (1/2)*heat_const_O2[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_O2[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_O2[2][0]*(T_f_in**4-T_f_out**4)

        cps_out = np.array([int_N2,int_CO2,int_H2O,int_O2])
        cp_f = np.dot(cps_out,mass_conc)/Mm_af #J/kg_fumée


        A_variables = np.array([heat_const_N2[2][3],heat_const_CO2[2][3],heat_const_H2O[2][3],heat_const_O2[2][3]])
        A_var = np.dot(A_variables,mass_conc)/Mm_af


        T_f_out_final = ((cp_f*Mflow_f_in - Q)/(A_var*Mflow_f_in))+T_f_in
        error = abs(T_f_out_final-T_f_out)
        T_f_out = T_f_out_final
        print("error",error)
        if iter==50 :
            print("Le heat exchanger ne converge peut-être pas")
            T_f_out = T_f_out_secours
    print("Nombres d'itération : ",iter)
    print("T_f_out : ",T_f_out-273.15,'C')


    cps_out = np.array([cpN2(T_f_out),cpCO2(T_f_out),cpH2O(T_f_out),cpO2(T_f_out)]) #premiere approx
    cp_f = np.dot(cps_out,mass_conc)/Mm_af #J/kg_fumée

    if Mflow_air_in*cp_air(T_air_in)>Mflow_f_in*cp_f : #je dois prendre quelle température pour cp_f et cp_air?
        Deltag_T = T_f_in - T_air_out
        Deltap_T = T_f_out - T_air_in
        DeltaM_T = (Deltag_T - Deltap_T)/np.log(Deltag_T/Deltap_T)
        S = Q/(U*DeltaM_T)
    if Mflow_air_in*cp_air(T_air_in)<Mflow_f_in*cp_f :
        Deltag_T = T_f_out - T_air_in
        Deltap_T = T_f_in - T_air_out
        DeltaM_T = (Deltag_T - Deltap_T)/np.log(Deltag_T/Deltap_T)
        S = Q/(U*DeltaM_T)

    print("S : ",S)

    """
    MEAN TEMPERATURE CALCULATIONS
    """

    cps_T = np.array([janaf_integrate(cpN2_T,T_f_out,T_f_in,dt),janaf_integrate(cpCO2_T,T_f_out,T_f_in,dt),janaf_integrate(cpH2O_T,T_f_out,T_f_in,dt),janaf_integrate(cpO2_T,T_f_out,T_f_in,dt)])
    cp_f_T = np.dot(cps_T,mass_conc)/Mm_af

    cpTout = np.array([janaf_integrate(cpN2,273.15,T_f_out,dt),janaf_integrate(cpCO2,273.15,T_f_out,dt),janaf_integrate(cpH2O,273.15,T_f_out,dt),janaf_integrate(cpO2,273.15,T_f_out,dt)])
    cpfTout = np.dot(cpTout,mass_conc)/Mm_af

    cpTin = np.array([janaf_integrate(cpN2,273.15,T_f_in,dt),janaf_integrate(cpCO2,273.15,T_f_in,dt),janaf_integrate(cpH2O,273.15,T_f_in,dt),janaf_integrate(cpO2,273.15,T_f_in,dt)])
    cpfTin = np.dot(cpTin,mass_conc)/Mm_a

    Tfm = (-cpfTout+cpfTin)/cp_f_T #[K] mean temperature des fumées
    #print(Tfm)

    ###################################

    Tam = (janaf_integrate(cp_air,273.15,T_air_out,dt) - janaf_integrate(cp_air,273.15,T_air_in,dt))/(janaf_integrate(cp_air_T,T_air_in,T_air_out,dt))
    #print(Tam) #[K] mean temperature air

    """
    En abscence de termes dissipatifs l'efficacité exergétique du transfert de chaleur vaut:
    """

    eta_transex = ((Tam - 288.15)*Tfm)/(Tam*(Tfm-288.15))
    print("eta_transex :",eta_transex)

    outputs = GT_comb_arg.exchanger_output();
    outputs.T_air_out = T_air_out-273.15 #[°C]
    outputs.T_f_out = T_f_out -273.15 #°C
    outputs.Q = Q/1000#[kJ]
    outputs.eta_transex = eta_transex
    outputs.Surf = S
    return outputs
sol = heatexchanger(GT_comb_arg.exchanger_input(U =3),480)
