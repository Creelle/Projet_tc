from matplotlib import pyplot as plt
import numpy as np

def janaf_integrate(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt)
def function(T):
    return 2*T

print(janaf_integrate(function,0,2,0.01)))
