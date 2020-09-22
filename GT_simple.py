
from thermochem import janaf
db = janaf.Janafdb();
import numpy;

import GT_arguments as GT_arg;
import combustionGT as comb;


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
    eta_PiC = arg_in.eta_PiC;
    if eta_PiC ==-1.:
        eta_PiC = 0.9;#max temperature = 1050°C
    eta_PiT = arg_in.eta_PiT;
    if eta_PiT ==-1.:
        eta_PiT = 0.9;#max temperature = 1050°C
    
    
    ## preliminary data (air)
    # ======================
    
    # Your job
    
    # cp air at 15°C (298K): [kJ/mol/K]
    O2 = db.getphasedata('O2','g');
    print(O2.cp(298))


   
    ## cycle definition
    # =================

    comb.combustionGT(GT_arg.comb_input())

    ## define output arguments
    # ======================
    outputs = GT_arg.GT_outputs();
    outputs.eta[1] = 0.35;
    
    # Your job
    
    return outputs;




#tests    
GT_simple_outputs = GT_simple(GT_arg.GT_input());
