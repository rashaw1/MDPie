###################################################
#                                                 #
#  MDPie by Robert A Shaw (2016)                  #
#                                                 #
###################################################
import numpy as np
import time
import ljlib

# Set cutoff distance and maxiter
#RHO_VAR = 0.125
#T_VAR = 2.0
rcut = 3.0
maxiter = 100000
equilib_iter = 10000 

# Initialise N atoms onto a cubic lattice, with boxes of length box
def initPositions(N, box):
    positions = np.zeros((N, 3), float)
    lattice = int(N**(1.0/3.0) + 1.0) # number of lattice sites
    r = box*(np.arange(lattice, dtype=float)/lattice - 0.5)
    # loop through lattice placing each atom
    i = 0
    for x in r:
        for y in r:
            for z in r:
                positions[i] = np.array([x, y, z], float)
                i += 1
                if i >= N:
                    return positions
    return positions

# Rescale the velocities to desired temperature, T
def rescale(velocities, T):
    # give zero momentum
    velocities = velocities - velocities.mean(axis=0)
    # Calculate kinetic energy
    ekin = 0.5*np.sum(velocities*velocities)
    scalefactor = np.sqrt(1.5*len(velocities)*T/ekin)
    velocities = velocities*scalefactor
    return velocities

# initialise velocities randomly
def initVelocities(N, T):
    velocities = np.random.rand(N, 3)
    velocities = rescale(velocities, T)
    return velocities

# test run
def md_run(RHO_VAR, T_VAR):
    N = 10000
    rho = RHO_VAR
    box = (N/rho)**(1.0/3.0)
    print 'Box size: ', box ,'\n'
    T = T_VAR
    dt = 0.01
    epot = 0.0
    ekin = 0.0
    freq = 20.0
    nthreads = 1
    if N > 100:
        nthreads = 2**(np.floor(np.log10(N))-1);
    print 'Using ', nthreads, ' threads\n'
    
        
    rescaleGap = 50
    displayRate = 0.1
    printGap = 20
    
    positions = initPositions(N, box)
    velocities = initVelocities(N, T)
    avg_pressure = 0.0
    avg_temp = 0.0
    rc3 = 1.0/(rcut**3)
    ptail = (16.0/3.0)*np.pi*((2.0/3.0)*(rc3**3) - rc3) 
    
    xyzfile = open('pymd_test.xyz', 'w')
    
    # Make the LJContainer
    ljobject = ljlib.PyLJContainer(box, rcut, dt)
    for i in range(N):
        x, y, z = positions[i][0], positions[i][1], positions[i][2]
        vx, vy, vz = velocities[i][0], velocities[i][1], velocities[i][2]
        ljobject.addParticle(x, y, z, vx, vy, vz)
    ljobject.forcesEnergies(nthreads)
    
    startTime = time.time()
    lastTime = time.time()
    i = 0
    while i < maxiter:
        # propagate positions/velocities/accelerations
        ljobject.integrate(nthreads)

        # rescale velocities if necessary
        if i % rescaleGap == 0:
            ljobject.andersen(T, freq)

        epot = ljobject.getEPot()
        ekin = ljobject.getEKin()
        virial = ljobject.getVirial()
        pressure = rho*(T - virial/float(3*N) ) + ptail*(rho**2)
        
        # write to screen if necessary
        if time.time() - lastTime > displayRate:
            print "%7d %11.4f %11.4f %11.4f %11.4f %11.4f" % (i, epot + ekin, epot, ekin, ekin/(1.5*N), pressure)
            lastTime = time.time()

        if i > equilib_iter:
            avg_pressure = avg_pressure + pressure
            avg_temp = avg_temp + ekin/(1.5*N)
            if i % printGap == 0:
                # Write xyz coordinates out
                xyzfile.write('%d\n XYZ\n' % (N))
                for j in range(N):
                    jpos = ljobject.getPos(j)
                    xyzfile.write('C %11.4f %11.4f %11.4f\n' % (jpos[0], jpos[1], jpos[2]))
                            
        i = i + 1
        
    epot = ljobject.getEPot()
    ekin = ljobject.getEKin()
    virial = ljobject.getVirial()
    pressure = rho*(T - virial/float(3*N) ) + ptail*(rho**2)
    print "%7d  %11.4f  %11.4f  %11.4f %11.4f %11.4f" % (i, epot+ekin, epot, ekin, ekin/(1.5*N), pressure)

    endTime = time.time()
    print "Total time: %.1f s" % (endTime - startTime)

    print "Density = %11.4f\nAverage pressure = %11.4f\nAverage temperature = %11.4f" % ( rho, avg_pressure/float(maxiter-2000), avg_temp/float(maxiter-2000))
    xyzfile.close()
    
if __name__ == '__main__':
    Tlist = [0.9]
    rholist = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    for TEMP in Tlist:
        for RHO in rholist:
            md_run(RHO, TEMP)


