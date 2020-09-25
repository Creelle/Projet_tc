# Install thermochem:
# pip install thermchem

from thermochem import combustion
from thermochem import janaf
from thermochem import burcat
db = janaf.Janafdb()
mixdb = burcat.Elementdb()
mix = mixdb.getmixturedata([("O2 REF ELEMENT", 20.9476), ("N2 REF ELEMENT",78.084), ("CO2", 0.0319), ("AR REF ELEMENT", 0.9365), ])
fuel=mixdb.getelementdata("CH4   RRHO")

combObject=combustion.SimpleCombustor(fuel,2,mixdb)
print(combObject.adiabatic_flame_temp(298.15+600))
O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
# CH4=db.getphasedata('CH4',phase ='g');
# print(CH4)

print(O2)

#heating value of di-oxygen at 298 K : 29.375 J/mol/K
print(O2.cp(298))

#heating value of di-nitrogen at 298 K : 29.375 J/mol/K
print(N2.cp(298))

#heating value of dioxygen at 298 K : 29.375 J/mol/K
print(CO2.cp(298))

#heating value of dioxygen at 298 K : 29.375 J/mol/K
print(H2O.cp(298))
print(H2O.cv(298))

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
    deltah = Cp_a*(T-298.15)
    return Cp_a,deltah;

#heating value of dioxygen at 298 K : 1.012 kJ/kg/K
print(air_mixture(298))
print(air_mixture(600))
