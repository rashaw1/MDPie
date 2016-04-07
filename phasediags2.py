import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.05, 0.35, 13)
y2 = [0.0938, 0.1365, 0.1776, 0.2168, 0.2558, 0.2915, 0.3275, 0.3732, 0.4029, 0.4477, 0.4849, 0.5318, 0.5796]
y09 = [0.0309, 0.0335, 0.0336, 0.0222, 0.0165, 0.005, -0.0073, -0.0143, -0.0337, -0.052, -0.0737, -0.0844, -0.1148]
y31 = [0.1542, 0.2317, 0.3089, 0.3897, 0.4716, 0.5551, 0.6417, 0.7322, 0.8387, 0.9361, 1.0554, 1.1745, 1.3204]
y265 = [0.1298, 0.1929, 0.2561, 0.3176, 0.3825, 0.448, 0.5195, 0.5857, 0.668, 0.7378, 0.8244, 0.9095, 1.0248]
y135 = [0.0576, 0.0797, 0.0972, 0.1124, 0.1228, 0.1291, 0.1362, 0.1433, 0.1426, 0.1507, 0.1471, 0.1379, 0.1422]
y125 = [0.0521, 0.0706, 0.085, 0.0932, 0.1018, 0.1054, 0.1045, 0.1072, 0.1044, 0.1026, 0.0925, 0.0881, 0.0793]

ylist = [y2, y09, y31, y265, y135, y125]
for y in ylist:
    plt.plot(x, y, 'k-')

plt.xlim([0.05, 0.35])
plt.ylim([-0.2, 1.2])
plt.xlabel(r'$\rho$')
plt.ylabel('P')
plt.yticks(np.arange(-0.2, 1.4, 0.2))
plt.show()

