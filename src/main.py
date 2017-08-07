

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
numberSimulationSteps = 1E5
numberElectronTraps   = 100
electronTravelRange   = 150E-9
rareEarthIndex        = [(i,0) for i in range(-50, 50+1)]
centerElectronTraps   = 0.0
spanElectronTraps     = 3.0E-6
pumpAmplitude		  = np.arange(0.01, 0.1, 0.01)
stedAmplitude		  = np.arange(5.0, 51.0, 1.0)
saturationAveraging   = 0.0001*numberSimulationSteps
					  # gamma, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE
crossSections 		  = [0.2,     1.0,       15.0,          2.0,       1.0]

rootPath = "C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/2D/"
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
	pass #raise ValueError('Path exists.')
else:
	os.makedirs(path)

simulatorList = list()
groundAverage  = list()
excitedAverage = list()
simResult = dict()

#--------------------------------------------------------------------------
# now simulate
#--------------------------------------------------------------------------

for sa in stedAmplitude:
	pa = pumpAmplitude
	groundAverage  = list()
	excitedAverage = list()
	simResult = dict()
	for val in rareEarthIndex:
		start_time = timeit.default_timer()
		s = SolidStateStedSimulator(nSimSteps = numberSimulationSteps,
									nET       = numberElectronTraps,
									eTR       = electronTravelRange,
									posRE     = val,
									centerET  = centerElectronTraps,
									spanET    = spanElectronTraps,
									pumpAmpl  = pa,
									stedAmpl  = sa,
									satAvg    = saturationAveraging,
									savePath  = path,
									cs        = crossSections)
		s.setupSimulation()

		s.start()
		s.join()
		#s.visualize()
		stop_time = timeit.default_timer()
		print "total runtime: %.1f s"%(stop_time - start_time), "pumpAmpl: %.2f, stedAmpl: %.2f"%(pa, sa)

		#path = "C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/"
		#pickle.dump(s, open("%sREidx_%d_pump_%.2f_sted_%.2f.pys"%(path, val, pa, sa), "wb"))

		groundAverage.append(np.average(np.array_split(s.electronSystems.groundStateEvolution._y, 2)[1]))
		excitedAverage.append(np.average(np.array_split(s.electronSystems.excitedStateEvolution._y, 2)[1]))

	# plot result (PSF)
	posStep = s.electronSystems.x[-1] - s.electronSystems.x[-2]
	xPos = np.array([i*posStep for i in range(len(rareEarthIndex))])
	#pp=s.pumpBeam.profile(s.electronTrapXPosition, s.electronTrapYPosition)
	plt.plot(xPos, np.array(excitedAverage)/np.max(excitedAverage))
	#plt.plot(xPos, pp/np.max(pp))
	plt.title("pump = %.2f, sted = %.2f"%(pa, sa))
	plt.savefig("%sPSF_pumpAmpl_%.2f_stedAmpl_%.2f.png"%(path, pa, sa))
	plt.close()

	# save simulation result (PSF)
	simResult['xPos'] = s.electronTrapXPosition
	simResult['yPos'] = s.electronTrapYPosition
	#simResult['pumpBeamProfile'] = pp
	simResult['excitedStateAverage'] = np.array(excitedAverage)
	simResult['pumpAmplitude'] = pa
	simResult['stedAmplitude'] = sa
	simResult['gamma'] = crossSections[0]
	simResult['sigPumpRE'] = crossSections[1]
	simResult['sigIonizeRE'] = crossSections[2]
	simResult['sigRepumpRE'] = crossSections[3]
	simResult['sigStedRE'] = crossSections[4]
	simResult['numberSimulationSteps'] = numberSimulationSteps
	simResult['numberElectronTraps'] = numberElectronTraps**2
	simResult['electronTravelRange'] = electronTravelRange
	pickle.dump(simResult, open("%sPSF_pumpAmpl_%.2f_stedAmpl_%.2f.pys"%(path, pa, sa),'wb'))

