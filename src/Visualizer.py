

import matplotlib.pyplot as plt
import numpy as np

#==============================================================================
class Visualizer(object):

	#--------------------------------------------------------------------------
	def __init__(self):
		pass

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	def visualize(self, pumpBeam, stedBeam, electronicSystems):

		electronicSystemXPositions = list()
		for es in electronicSystems:
			electronicSystemXPositions.append(es.x)

		electronicSystemXPositions = np.array(electronicSystemXPositions)

		a=pumpBeam.profile(electronicSystemXPositions)
		b=stedBeam.profile(electronicSystemXPositions)
		plt.plot(electronicSystemXPositions/1.0E-6, a)
		plt.plot(electronicSystemXPositions/1.0E-6, b)

		for es in electronicSystems:
			plt.plot(es.x/1.0E-6, 0.0009, 'ro')
			if es.isPopulated:
				plt.plot(es.x/1.0E-6, 0.001, 'g^')

		plt.show()
		

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
