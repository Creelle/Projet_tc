from matplotlib import pyplot as plt
import numpy as np

x= np.linspace(0,10,1000)
plt.rcParams.update({'font.size': 24})


fig, ax= plt.subplots()

ax.plot(x,x,'-b')
ax.set_xlabel("Number of integrals done")
ax.set_ylabel("Intuition")
plt.grid(True)
plt.xticks([])
plt.yticks([])
plt.show()
