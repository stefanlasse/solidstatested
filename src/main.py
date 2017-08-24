

import timeit
import pickle
import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['savefig.directory'] = os.chdir(os.path.dirname(__file__))
mpl.rcParams['font.size'] = 16

from Simulator import SolidStateStedSimulator

#--------------------------------------------------------------------------
# configuration part
#--------------------------------------------------------------------------
numberSimulationSteps = 5E5
numberElectronTraps   = 50

rareEarthXposition    = np.linspace(-2.5E-7, 2.5E-7, 61)
rareEarthYposition    = np.zeros(rareEarthXposition.size)

electronTrapXposition = np.linspace(-5E-7, 5E-7, numberElectronTraps)
electronTrapYposition = np.linspace(-5E-7, 5E-7, numberElectronTraps)

pumpAmplitude		  = np.linspace(0.05, 0.1, 6)
stedAmplitude		  = np.linspace(1.0, 20.0, 20)

					  #  gammaRE, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE
crossSections 		  = [  0.2,      2.0,       7.0,         5.0,        3.0]

electronTravelRange   = 101E-9

rootPath = "D:/STED_sim/2D/"
path = "%sgamma_%.2f_sigPumpRE_%.2f_sigIonizeRE_%.2f_sigRepumpRE_%.2f_sigStedRE_%.2f/"%(rootPath,			\
																						crossSections[0],	\
																						crossSections[1],	\
																						crossSections[2],	\
																						crossSections[3],	\
																						crossSections[4])



#--------------------------------------------------------------------------
# some internals
#--------------------------------------------------------------------------
rareEarthCoordinates = np.vstack((rareEarthXposition, rareEarthYposition)).T
electronTrapCoordinates = np.array([[x,y] for x in electronTrapXposition for y in electronTrapYposition])

if os.path.exists(path):
	raise ValueError('Path exists. Simulation already done.')
else:
	os.makedirs(path)

#--------------------------------------------------------------------------
# now simulate
#--------------------------------------------------------------------------

for pa in pumpAmplitude:
	for sa in stedAmplitude:

		groundAverage  = list()
		excitedAverage = list()
		simResult = dict()

		for rePos in rareEarthCoordinates:
			start_time = timeit.default_timer()

			s = SolidStateStedSimulator(nSimSteps = numberSimulationSteps, savePath = path)

			s.setupSimulation(REx=rePos[0], REy=rePos[1],
							  ETx=electronTrapCoordinates[:,0], ETy=electronTrapCoordinates[:,1],
							  pumpAmpl=pa, stedAmpl=sa,
							  cs=crossSections, eTR=electronTravelRange)
			s.start()
			s.join()

			stop_time = timeit.default_timer()
			print "total runtime: %.1f s\n"%(stop_time - start_time), "pumpAmpl: %.2f, stedAmpl: %.2f\n"%(pa, sa), "REpos: [%.3g, %.3g] \n"%(rePos[0], rePos[1])

			groundAverage.append(np.average(np.array_split(s.evolutionRecoders[0]._g, 2)[1]))
			excitedAverage.append(np.average(np.array_split(s.evolutionRecoders[0]._e, 2)[1]))

		# plot result (PSF)
		plt.plot(rareEarthXposition, np.array(excitedAverage)/np.max(excitedAverage))
		plt.title("pump = %.2f, sted = %.2f"%(pa, sa))
		#plt.show()
		plt.savefig("%sPSF_pumpAmpl_%.2f_stedAmpl_%.2f.png"%(path, pa, sa))
		plt.close()

		# save simulation result (PSF)
		simResult['RExPos'] = rareEarthXposition
		simResult['REyPos'] = rareEarthYposition
		simResult['groundStateAverage'] = np.array(groundAverage)
		simResult['excitedStateAverage'] = np.array(excitedAverage)
		simResult['pumpAmplitude'] = pa
		simResult['stedAmplitude'] = sa
		simResult['gammaRE'] = crossSections[0]
		simResult['sigPumpRE'] = crossSections[1]
		simResult['sigIonizeRE'] = crossSections[2]
		simResult['sigRepumpRE'] = crossSections[3]
		simResult['sigStedRE'] = crossSections[4]
		simResult['numberSimulationSteps'] = numberSimulationSteps
		simResult['numberElectronTraps'] = numberElectronTraps**2
		simResult['electronTravelRange'] = electronTravelRange
		pickle.dump(simResult, open("%sPSF_pumpAmpl_%.2f_stedAmpl_%.2f.pys"%(path, pa, sa),'wb'))

