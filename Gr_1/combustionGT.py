## Project LMECA2150-Thermal cycle
# Material related to the combustion of the Gas turbine
#
# Author: Paolo Thiran & Gauthier Limpens
# Version: 2020
#
# This script can be modified entirely by the students.


from thermochem import janaf
db = janaf.Janafdb();
import numpy;

import GT_arguments as GT_arg;
import GTcomb_arguments as GTcomb_arg;


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
       -options.r     [-] : Comperssion ratio
                            chamber
       -options.T_3   [°C] : Temperature after combustion (before turbine)
       -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for compression
       -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for expansion
    """
    
    # Your job
    #print(comb_input.lambda_comb);
    
    ## define output arguments
    # ======================
    outputs = GTcomb_arg.comb_output();
    outputs.T_out = 1543;
    
    # Your job
    
    return outputs;

### test
## combustionGT(GT_arg.comb_input())
