from matplotlib import pyplot as plt
import numpy as np

def janaf_integrate(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt)
def function(T):
    return 2*T


print(np.arange(100,0,10))
