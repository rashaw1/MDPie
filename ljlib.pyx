# distutils: language = c++
# distutils: sources = ljforces.cpp

from libcpp.vector cimport vector

cdef extern from "ljforces.hpp" namespace "lj":
	cdef cppclass LJContainer:
		LJContainer(double, double, double) except +
		int N
		vector[double] positions, velocities, forces
		double epot, ekin, virial, box, rcutoff, dt
		int getN()
		double getEPot()
		double getEKin()
		double getVirial()
		vector[double] getPos(int)
		vector[double] getVel(int)
		vector[double] getForces(int)
		void setPos(int, double, double, double)
		void setVel(int, double, double, double)
		void setConsts(double, double, double)
		void addParticle(double, double, double, double, double, double)
		void removeParticle()
		void forcesEnergies(int)
		void forcesThread(int, int, vector[vector[double]], double, double, vector[vector[double]], double, int, double)
		void integrate(int)
		void andersen(double, double)
		
cdef class PyLJContainer:
	cdef LJContainer *thisptr
	def __cinit__(self, box, rcut, dt):
		self.thisptr = new LJContainer(box, rcut, dt)
		
	def getN(self):
		return self.thisptr.getN()
	def getEPot(self):
		return self.thisptr.getEPot()
	def getEKin(self):
		return self.thisptr.getEKin()
	def getVirial(self):
		return self.thisptr.getVirial()
		
	def getPos(self, int i):
		return self.thisptr.getPos(i)
	def getVel(self, int i):
		return self.thisptr.getVel(i)
	def getForces(self, int i):
		return self.thisptr.getForces(i)
		
	def setPos(self, int i, x, y, z):
		self.thisptr.setPos(i, x, y, z)
	def setVel(self, int i, x, y, z):
		self.thisptr.setVel(i, x, y, z)
	def setConsts(self, box, rcut, dt):
		self.thisptr.setConsts(box, rcut, dt)
	
	def addParticle(self, x, y, z, vx, vy, vz):
		self.thisptr.addParticle(x, y, z, vx, vy, vz)
	def removeParticle(self):
		self.thisptr.removeParticle()
	
	def forcesEnergies(self, int nthreads):
		self.thisptr.forcesEnergies(nthreads)
	def integrate(self, int nthreads):
		self.thisptr.integrate(nthreads)
	def andersen(self, T, freq):
		self.thisptr.andersen(T, freq)
