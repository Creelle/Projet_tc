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
     OPTIONS is a structure containing :
       -options.k_mec [-] : Shaft losses
       -options.T_0   [°C] : Reference temperature
       -options.T_ext [°C] : External temperature
       -options.r     [-] : Compression ratio
       -options.k_cc  [-] : Coefficient of pressure losses due to combustion
                            chamber
       -options.T_3   [°C] : Temperature after combustion (before turbine)
       -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for compression
       -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for expansion
     DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
              do not plot."""
    def __init__(self, Pe = 50e3,#50 MW
                     k_mec = 0,
                     T_0 = 15.,#°C
                     T_ext =288.15,
                     r =10.,
                     k_cc =-1.,
                     T3 =1050,
                     eta_PiC =0.9,
                     eta_PiT =0.9,
                     Display =-1.):
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
       -eta(3) : eta_cyclex, cycle exergy efficiency
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
      "FIG(1) = figure;
      plot(x,y1);
      [...]
       FIG(2) = figure;
      plot(x,y2);
      [...]"
      Your vector FIG will contain all the figure plot during the run of this
      code (whatever the size of FIG).
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
                     h_in = 650,# enthalpy[kJ/kg_air]
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
