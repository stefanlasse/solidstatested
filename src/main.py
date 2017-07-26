

import timeit

from Simulator import SolidStateStedSimulator


tl = list()

for pa in [0.005, 0.01, 0.015, 0.02, 0.025]:
	s = SolidStateStedSimulator(nSimSteps=5E6, nET=300, eTR=25E-9, posRE=0.0, centerET=0.0, spanET=2E-6,
								pumpAmpl=pa, stedAmpl=0.0)
	tl.append(s)


start_time = timeit.default_timer()
for s in tl:
	s.setupSimulation()
	s.start()


for s in tl:
	s.join()


stop_time = timeit.default_timer()
print "total runtime: %.3f s"%(stop_time - start_time)

#s.visualize()