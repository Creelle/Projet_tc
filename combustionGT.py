
from thermochem import janaf
from thermochem import combustion
from thermochem import burcat
db = janaf.Janafdb();
burcatdb = burcat.Elementdb()
# mix = burcatdb.getmixturedata([("O2 REF ELEMENT", 20.9476), ("N2 REF ELEMENT",78.084), ("CO2", 0.0319), ("AR REF ELEMENT", 0.9365), ])

import numpy as np;

import GT_arguments as GT_arg;

def combustionGT(comb_input):
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
       -options.r     [-] : Compression ratio
                            chamber
       -options.T_3   [°C] : Temperature after combustion (before turbine)
       -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for compression
       -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for expansion
    """

    lambda_comb = comb_input.lambda_comb
    x_O2a = comb_input.x_O2a
    x_N2a = comb_input.x_N2a
    T_in  = comb_input.T_in
    h_in  = comb_input.h_in
    LHV   =comb_input.LHV
    # CH4 + 2 *lambda * (O2 + coeff*N2) <=> CO2+2*H2O+ 2*lambda*coeff*N2 + (2*lambda-2)*O2

    fuel=burcatdb.getelementdata("CH4   RRHO")  # choose the right fuel
    combObject=combustion.SimpleCombustor(fuel,1/lambda_comb,burcatdb)

    # donne une liste avec les coeff stochio pour les produits de la combustion dans cette ordre N2 - CO2 - H2O - O2
    coeff_stochio = combustion.balance(fuel,1,1/lambda_comb)[1]
    print(combustion.balance(fuel,1,1/lambda_comb))
    list_coeff = np.asarray(list(coeff_stochio.values()))

    total_n = sum(list_coeff) # nombre de moles total
    molar_conc = list_coeff/total_n # concentration des elements mol_co2/mol_t
    molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
    Mm_af = sum(molar_conc*molar_mass) #somme ponderé des masse molaire

    #  OUTPUTS
    #  R_f   : ideal gas constant (R*) for specific gas (R/Mm_f) [kJ/kg/K]
    #  m_O2f : mass fraction of oxygen in exhaust gases [-]
    #  m_N2f : mass fraction of nitrogen in exhaust gases [-]
    #  m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
    #  m_H2Of : mass fraction of water steam in exhaust gases [-]
    #  T_out  : outlet gas temperature [K]

    outputs = GT_arg.comb_output();
    outputs.T_out = combObject.adiabatic_flame_temp(273.15+T_in); #[K]
    outputs.R_f = 8.31/1000/Mm_af # [kJ/kg/K]
    outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = molar_conc*molar_mass/Mm_af #[-]

    return outputs;

### test

sol =combustionGT(GT_arg.comb_input(lambda_comb = 2))
print(sol.T_out)
