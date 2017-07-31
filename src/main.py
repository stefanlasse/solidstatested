

import timeit
from Simulator import SolidStateStedSimulator
import pickle

numberSimulationSteps = 1E5
numberElectronTraps   = 300
electronTravelRange   = 41E-9
rareEarthIndex        = 0
centerElectronTraps   = 0.0
spanElectronTraps     = 3.0E-6

simulatorList = list()

s = SolidStateStedSimulator(nSimSteps=numberSimulationSteps,
							nET=numberElectronTraps,
							eTR=electronTravelRange,
							posRE=rareEarthIndex,
							centerET=centerElectronTraps,
							spanET=spanElectronTraps)
s.setupSimulation()

start_time = timeit.default_timer()
s.start()
s.join()
stop_time = timeit.default_timer()
print "total runtime: %.1f s"%(stop_time - start_time)

s.visualize()