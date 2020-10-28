## Project LMECA2150-Thermal cycle
# Material related to the arguments related to the gas turbine combustion
#
# Author: Paolo Thiran & Gauthier Limpens
# Version: 2020
#
# This script can be modified entirely by the students.

import numpy;


class comb_input:
    """
     combustion(lambda,x_O2a,x_N2a) computes the combustion of methane with air given the excess air (lambda).
     It returns the molar fraction of the different 
    
     INPUTS 
     lambda_comb : excess air
     x_O2a : molar fraction of oxygen concentration in air
     x_N2a : molar fraction of nitrogen concentration in air
     T_in  : gas temperature at the inlet
     h_in  : enthalpy of the gas at the inlet
     LHV  : Fuel 'Low Heating Value'. CH4 here [kJ/kg_CH4] 
    
    """
    def __init__(self, lambda_comb = 2,#excess air
                     x_O2a = 0.21,# molar fraction 
                     x_N2a = 0.79,# molar fraction 
                     T_in = 600,#°C
                     h_in = 650,# enthalpy
                     LHV =50000):# [kJ/kg_CH4] 
        self.lambda_comb = lambda_comb;
        self.x_O2a = x_O2a;
        self.x_N2a = x_N2a ;
        self.T_in = T_in;
        self.h_in = h_in;
        self.LHV = LHV;
        
class comb_output:
    """
     combustion(lambda,x_O2a,x_N2a) computes the combustion of methane with air given the excess air (lambda).
     It returns the molar fraction of the different 
    
     OUTPUTS 
     R_f   : ideal gas constant (R*) for specific gas (R/Mm_f) [kJ/kg/K]
     m_O2f : mass fraction of oxygen in exhaust gases [-]
     m_N2f : mass fraction of nitrogen in exhaust gases [-]
     m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
     m_H2Of : mass fraction of water steam in exhaust gases [-]
     T_out  : outlet gas temperature [K]
    
    """
    def __init__(self, R_f = -1.,#50 MW
                     m_O2f = 0.1,# molar fraction 
                     m_N2f = 0.4,# molar fraction 
                     m_CO2f = 0.2,#°C
                     m_H2Of = 0.3,# enthalpy
                     T_out = -1.):# [kJ/kg_CH4] 
        self.R_f = R_f;
        self.m_O2f = m_O2f;
        self.m_N2f = m_N2f ;
        self.m_CO2f = m_CO2f;
        self.m_H2Of = m_H2Of;
        self.T_out = T_out;
        
        
