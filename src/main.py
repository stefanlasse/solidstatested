

import timeit
import numpy as np
import os
from multiprocessing import freeze_support

from PointSpreadFunction import PointSpreadFunction

#--------------------------------------------------------------------------
# configuration part
#--------------------------------------------------------------------------
numberSimulationSteps = 1E5
numberElectronTraps   = 50

rareEarthXposition    = np.linspace(-2.5E-7, 2.5E-7, 3)
rareEarthYposition    = np.zeros(rareEarthXposition.size)

electronTrapXposition = np.linspace(-5E-7, 5E-7, numberElectronTraps)
electronTrapYposition = np.linspace(-5E-7, 5E-7, numberElectronTraps)

pumpAmplitude		  = np.linspace(0.05, 0.1, 1)
stedAmplitude		  = np.linspace(1.0, 20.0, 1)

					  #  gammaRE, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE
crossSections 		  = [  0.2,      2.0,       7.0,         5.0,        3.0]

electronTravelRange   = 101E-9

rootPath = "C:/Users/lasse/Documents/projects/solidstatested/src/"
path = "%sgamma_%.2f_sigPumpRE_%.2f_sigIonizeRE_%.2f_sigRepumpRE_%.2f_sigStedRE_%.2f/"%(rootPath,			\
																						crossSections[0],	\
																						crossSections[1],	\
																						crossSections[2],	\
																						crossSections[3],	\
																						crossSections[4])

#--------------------------------------------------------------------------
# some internals
#--------------------------------------------------------------------------
if os.path.exists(path):
	pass #raise ValueError('Path exists. Simulation already done.')
else:
	os.makedirs(path)


rareEarthCoordinates = np.vstack((rareEarthXposition, rareEarthYposition)).T
electronTrapCoordinates = np.array([[x,y] for x in electronTrapXposition for y in electronTrapYposition])

#--------------------------------------------------------------------------
# now simulate
#--------------------------------------------------------------------------

if __name__ == '__main__':
	freeze_support()

	for pa in pumpAmplitude:
		for sa in stedAmplitude:

			start_time = timeit.default_timer()

			psf = PointSpreadFunction(numberSimulationSteps, rareEarthCoordinates, electronTrapCoordinates, pa, sa, crossSections, electronTravelRange, path)
			psf.start()
			psf.join()

			stop_time = timeit.default_timer()
			print "total runtime: %.1f s\n"%(stop_time - start_time)
			print "pump=%.2f, sted=%.1f"%(pa, sa)
			print ""
