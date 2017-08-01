

import timeit
from Simulator import SolidStateStedSimulator
import pickle
import numpy as np

import matplotlib.pyplot as plt

numberSimulationSteps = 1E5
numberElectronTraps   = 300
electronTravelRange   = 41E-9
rareEarthIndex        = [i for i in range(-100, 100+1)]
centerElectronTraps   = 0.0
spanElectronTraps     = 3.0E-6
pumpAmplitude		  = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 50.])
stedAmplitude		  = np.array([0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 50.])

# note for STED: pumpAmpl = 0.05
#                stedAmpl = np.array([0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 50.])


simulatorList = list()

groundAverage  = list()
excitedAverage = list()
simResult = dict()




pumpAmplitude		  = np.array([0.08, 0.13, 0.17, 0.3, 0.4, 7., 15., 30., 40., 70., 100.])
stedAmplitude		  = 0.0
#sweep pump amplitude
for pa in pumpAmplitude:
	sa = stedAmplitude
	groundAverage  = list()
	excitedAverage = list()
	simResult = dict()
	for val in rareEarthIndex:
		start_time = timeit.default_timer()
		s = SolidStateStedSimulator(nSimSteps=numberSimulationSteps,
									nET=numberElectronTraps,
									eTR=electronTravelRange,
									posRE=val,
									centerET=centerElectronTraps,
									spanET=spanElectronTraps,
									pumpAmpl=pa,
									stedAmpl=sa)
		s.setupSimulation()

		s.start()
		s.join()
		stop_time = timeit.default_timer()
		print "total runtime: %.1f s"%(stop_time - start_time), "pumpAmpl: %.2f, stedAmpl: %.2f"%(pa, sa)

		groundAverage.append(np.average(s.electronSystems.groundStateEvolution._y))
		excitedAverage.append(np.average(s.electronSystems.excitedStateEvolution._y))

	# plot result
	posStep = s.electronSystems.x[-1] - s.electronSystems.x[-2]
	xPos = np.array([i*posStep for i in range(len(rareEarthIndex))])
	pp=s.pumpBeam.profile(xPos-1E-6)
	plt.plot(xPos, np.array(excitedAverage)/np.max(excitedAverage))
	plt.plot(xPos, pp/np.max(pp))
	plt.title("pump = %.2f, sted = %.2f"%(pa, sa))
	plt.savefig("C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/PSF_pumpAmpl_%.2f_stedAmpl_%.2f.png"%(pa, sa))
	plt.close()

	# save simulation
	simResult['xPos'] = xPos
	simResult['pumpBeamProfile'] = pp
	simResult['excitedStateAverage'] = np.array(excitedAverage)
	simResult['pumpAmplitude'] = pa
	simResult['stedAmplitude'] = sa
	pickle.dump(simResult, open("C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/PSF_pumpAmpl_%.2f_stedAmpl_%.2f.pys"%(pa, sa), 'wb'))









rareEarthIndex        = [i for i in range(-40, 40+1)]
pumpAmplitude		  = 0.05
stedAmplitude		  = np.arange(4.0, 40.0, 4.0)
# sweep sted amplitude
for sa in stedAmplitude:
	pa = pumpAmplitude
	groundAverage  = list()
	excitedAverage = list()
	simResult = dict()
	for val in rareEarthIndex:
		start_time = timeit.default_timer()
		s = SolidStateStedSimulator(nSimSteps=numberSimulationSteps,
									nET=numberElectronTraps,
									eTR=electronTravelRange,
									posRE=val,
									centerET=centerElectronTraps,
									spanET=spanElectronTraps,
									pumpAmpl=pa,
									stedAmpl=sa)
		s.setupSimulation()

		s.start()
		s.join()
		stop_time = timeit.default_timer()
		print "total runtime: %.1f s"%(stop_time - start_time), "pumpAmpl: %.2f, stedAmpl: %.2f"%(pa, sa)

		groundAverage.append(np.average(s.electronSystems.groundStateEvolution._y))
		excitedAverage.append(np.average(s.electronSystems.excitedStateEvolution._y))

	# plot result
	posStep = s.electronSystems.x[-1] - s.electronSystems.x[-2]
	xPos = np.array([i*posStep for i in range(len(rareEarthIndex))])
	pp=s.pumpBeam.profile(xPos-1E-6)
	plt.plot(xPos, np.array(excitedAverage)/np.max(excitedAverage))
	plt.plot(xPos, pp/np.max(pp))
	plt.title("pump = %.2f, sted = %.2f"%(pa, sa))
	plt.savefig("C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/PSF_pumpAmpl_%.2f_stedAmpl_%.2f.png"%(pa, sa))
	plt.close()

	# save simulation
	simResult['xPos'] = xPos
	simResult['pumpBeamProfile'] = pp
	simResult['excitedStateAverage'] = np.array(excitedAverage)
	simResult['pumpAmplitude'] = pa
	simResult['stedAmplitude'] = sa
	pickle.dump(simResult, open("C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/PSF_pumpAmpl_%.2f_stedAmpl_%.2f.pys"%(pa, sa),'wb'))




