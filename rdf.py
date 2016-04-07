# Takes output of MD simulation and computes radial distribution function
# and plots a graph of it

import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d

FILENAME = "pymd_test.xyz"
OUTFILE = "ljfluid05.rdf"
NBOXES = 50 # Number of histogram boxes
BOXSIZE = 11.856311015 
RMAX = 5.0 # Maximum radius

# Open the file
infile = open(FILENAME, 'r')

# First line needs to be the number of particles
N = int(infile.readline())

# Set up histogram boxes
dr = RMAX/float(NBOXES)
hist = np.zeros(NBOXES)

# Loop through file
counter = 1
configs = 1
firstpos = np.zeros(3)
currpos = np.zeros(3)
for line in infile:
    if counter == N+2:
        counter = 0
        configs += 1
    else:
        parts = line.split()
        if counter == 2:
            postemp = [ float(x) for x in parts[1:4] ]
            firstpos[:] = postemp[:]
        elif counter > 2:
            postemp = [ float(x) for x in parts[1:4] ]
            currpos[:] = postemp[:]
            rij = currpos - firstpos
            r = np.sqrt(np.dot(rij, rij))
            i = 0
            while i < NBOXES: 
                if r < i*dr:
                    hist[i] += 1
                    i = N
                i += 1
    counter += 1

infile.close()

# Form g(r)
constant = 4.0*np.pi*(N/(float(BOXSIZE)**3))*configs/3.0

outfile = open(OUTFILE, 'w')

gr = np.zeros(NBOXES)
ri = np.zeros(NBOXES)
for i in range(1, NBOXES):
    ri[i] = i*dr
    shellsize = ((i+1)*dr)**3 - (i*dr)**3
    gr[i] = hist[i]/(constant*shellsize)
    outfile.write("%10.4f %10.4f\n" % (ri[i], gr[i]))
outfile.close()
    
# Plot

f = interp1d(ri, gr, kind='cubic')
xnew = np.linspace(0, max(ri), num=4*NBOXES, endpoint = True)
plt.plot(xnew, f(xnew), color='r')
plt.ylabel('g(r)')
plt.xlabel('r')
plt.show()
        
        
