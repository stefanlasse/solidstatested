

import timeit

from Simulator import SolidStateStedSimulator

s = SolidStateStedSimulator(nSimSteps=1E6, nET=300, eTR=41E-9, posRE=0, centerET=0.0, spanET=3E-6)

s.setupSimulation()

start_time = timeit.default_timer()
s.start()
s.join()
stop_time = timeit.default_timer()
print "total runtime: %.1f s"%(stop_time - start_time)

s.visualize()