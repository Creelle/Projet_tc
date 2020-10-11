# Install thermochem:
# pip install thermchem

from thermochem import combustion
from thermochem import janaf
from thermochem import burcat
db = janaf.Janafdb()

T0=288.15
O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
CH4=db.getphasedata('CH4',phase ='g');
C2H4 = db.getphasedata('C2H4',phase ='g');
print(O2)
print(N2)
print(CO2)
print(H2O)
# print(O2.hef([273.15,288.15,500]))
# print(N2.hef([273.15,288.15,500]))# kJ/mol
# print(CH4.hef([273.15,288.15,500]))
# #heating value of di-oxygen at 298 K : 29.375 J/mol/K
# print(O2.cp(298))



# Chemical properties:
Mm_O2 = 0.032;#kg/mol
Mm_N2 = 0.028;#kg/mol
conc_O2 = 0.21;# 21% in molar
conc_N2 = 0.79;# 79% in molar

def air_mixture(T):#kJ/kg/K
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp_a = cp_a/Mm_a/1000;#kJ/kg/K
    return Cp_a;

#heating value of dioxygen at 298 K : 1.012 kJ/kg/K
# print(air_mixture(298))
# print(air_mixture(600))
"""
exergie test
"""
print("dfsfsrgrg exergie test")
print(N2.hef(T0))
T = 400
e = (N2.hef(T)-N2.hef(T0))*1000-T0*(N2.S(T)-N2.S(T0))
