import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

NBOXES = 50

lj001 = np.loadtxt('ljfluid001.rdf')
lj005 = np.loadtxt('ljfluid005.rdf')
lj05 = np.loadtxt('ljfluid05.rdf')


f001 = interp1d(lj001[:,0], lj001[:,1], kind='cubic')
f005 = interp1d(lj005[:,0], lj005[:,1], kind='cubic')
f05 = interp1d(lj05[:,0], lj05[:,1], kind='cubic')
x001 = np.linspace(0.1, max(lj001[:,0]), num=4*NBOXES, endpoint = True)
x005 = np.linspace(0.1, max(lj005[:,0]), num=4*NBOXES, endpoint = True)
x05 = np.linspace(0.1, max(lj05[:,0]), num=4*NBOXES, endpoint = True)
plt.plot(x001, f001(x001), color='r', label='T=0.01')
plt.plot(x005, f005(x005), color='g', linestyle='dashed', label='T=0.05')
plt.plot(x05, f05(x05), color='b', linestyle='dotted', label='T=0.5')
plt.ylabel('g(r)')
plt.xlabel('r')
plt.legend(loc='upper right')
plt.show()
