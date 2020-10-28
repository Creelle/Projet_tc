## Project LMECA2150-Thermal cycle
# Material related to the combustion of the Gas turbine
#
# Author: Paolo Thiran & Gauthier Limpens
# Version: 2020
#
# !!! This script CANNOT be modified by the students.

import numpy;

class GT_input:
    """ GT Gas turbine modelisation
     GT(P_e,options,display) compute the thermodynamics states for a Gas
     turbine based on several inputs (given in OPTION) and based on a given 
     electricity production P_e. It returns the main results. It can as well
     plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
    
     INPUTS (some inputs can be dependent on others => only one of these 2 can
             be activated) Refer to Fig 3.1 from reference book (in english)
     P_E = electrical power output target [kW]
     k_mec [-] : Shaft losses 
     T_0   [°C] : Reference temperature
     T_ext [°C] : External temperature
     r     [-] : Compression ratio
     k_cc  [-] : Coefficient of pressure losses due to combustion
                 chamber
     T_3   [°C] : Temperature after combustion (before turbine)
     eta_PiC[-] : Intern polytropic efficiency (Rendement
                  polytropique interne) for compression
     eta_PiT[-] : Intern polytropic efficiency (Rendement
                  polytropique interne) for expansion
     DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
              do not plot."""
    def __init__(self, Pe = 50e3,#50 MW
                     k_mec = -1.,
                     T_0 = -1.,#°C
                     T_ext =10.,#°C
                     r =10.,
                     k_cc =-1.,
                     T3 =1050,
                     eta_PiC =0.9,
                     eta_PiT =0.9,
                     Display =0):
        self.Pe = Pe;
        self.k_mec = k_mec;
        self.T_0 = T_0 ;
        self.T_ext = T_ext;
        self.r = r;
        self.k_cc = k_cc;
        self.T3 = T3;
        self.eta_PiC = eta_PiC;
        self.eta_PiT = eta_PiT;
        self.Display = Display;
    
class GT_outputs:
    """
     OUPUTS : 
     ETA is a vector with :
       -eta(1) : eta_cyclen, cycle energy efficiency
       -eta(2) : eta_toten, overall energy efficiency
       -eta(3) : eta_cyclex, cycle exegy efficiency
       -eta(4) : eta_totex, overall exergie efficiency
       -eta(5) : eta_rotex, compressor-turbine exergy efficiency
       -eta(6) : eta_combex, Combustion exergy efficiency
       FYI : eta(i) \in [0;1] [-]
     DATEN is a vector with : 
       -daten(1) : perte_mec [kW]
       -daten(2) : perte_ech [kW]
     DATEX is a vector with :
       -datex(1) : perte_mec [kW]
       -datex(2) : perte_rotex [kW]
       -datex(3) : perte_combex [kW]
       -datex(4) : perte_echex  [kW]
     DAT is a matrix containing :
     dat = {T_1       , T_2       , T_3       , T_4; [°C]
            p_1       , p_2       , p_3       , p_4; [bar]
            h_1       , h_2       , h_3       , h_4; [kJ/kg]
            s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
            e_1       , e_2       , e_3       , e_4;};[kJ/kg]
     MASSFLOW is a vector containing : 
       -massflow(1) = m_a, air massflow [kg/s]
       -massflow(2) = m_c, combustible massflow [kg/s] 
       -massflow(3) = m_f, exhaust gas massflow [kg/s]
     
     COMBUSTION is a structure with :
       -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
       -combustion.e_c    : the combuistible exergie         [kJ/kg]
       -combustion.Lambda : the air excess                   [-]
       -combustion.Cp_g   : heat capacity of exhaust gas at 400 K [kJ/kg/K]
       -combustion.fum  : is a vector of the exhaust gas composition :
           -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
           -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
           -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
           -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 

     FIG is a vector of all the figure you plot. Before each figure, define a
     figure environment such as:  
      "fig1 = plt.figure(1);
      plt.plot(t, y1);
      [...]
       fig1 = plt.figure(1);
      plt.plot(t, y1);
      [...]
      fig=[fig1,fig2]"
      Your vector FIG will contain all the figure plot during the run of this
      code.
    """
    class Combustion:
        def __init__(self):
            self.LHV = 0.;
            self.e_c = 0.;
            self.Lambda = 0.;
            self.Cp_g = 0.;
            self.fum = numpy.zeros(4);
            
    def __init__(self, comb = Combustion()):
        self.eta = numpy.zeros(6);
        self.daten = numpy.zeros(2);
        self.datex = numpy.zeros(4);
        self.dat = numpy.zeros((5, 4));
        self.massflow = numpy.zeros(3);
        self.combustion = comb;
        self.fig = list();

