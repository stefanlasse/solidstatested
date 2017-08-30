

import timeit
import numpy as np
import os
from multiprocessing import freeze_support

from PointSpreadFunction import PointSpreadFunction

#--------------------------------------------------------------------------
# configuration part
#--------------------------------------------------------------------------
numberSimulationSteps = 5E5
numberElectronTraps   = 50

laserXposition = np.linspace(-2.5E-7, 2.5E-7, 63)
laserYposition = np.zeros(laserXposition.size)

rareEarthCoordinates = np.array([0.0, 0.0])

electronTrapXposition = np.linspace(-5E-7, 5E-7, numberElectronTraps)
electronTrapYposition = np.linspace(-5E-7, 5E-7, numberElectronTraps)

pumpAmplitude		  = np.linspace(0.05, 0.1, 1)
stedAmplitude		  = np.linspace(1.0, 20.0, 10)

					  #  gammaRE, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE
crossSections 		  = [  0.2,      2.0,       10.0,         5.0,        1.0]

electronTravelRange   = 101E-9

rootPath = "D:/STED_sim/test"
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


laserCoordinates = np.vstack((laserXposition, laserYposition)).T
electronTrapCoordinates = np.array([[x,y] for x in electronTrapXposition for y in electronTrapYposition])

#--------------------------------------------------------------------------
# now simulate
#--------------------------------------------------------------------------

if __name__ == '__main__':
	freeze_support()

	for pa in pumpAmplitude:
		for sa in stedAmplitude:

			start_time = timeit.default_timer()

			psf = PointSpreadFunction(numberSimulationSteps, rareEarthCoordinates, electronTrapCoordinates, pa, sa, laserCoordinates, crossSections, electronTravelRange, path)
			psf.start()
			psf.join()

			stop_time = timeit.default_timer()
			print "total runtime: %.1f s"%(stop_time - start_time)
			print "pump=%.2f, sted=%.1f"%(pa, sa)
			print ""
