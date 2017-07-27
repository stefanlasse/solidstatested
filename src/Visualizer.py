

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

		esXPos = electronicSystems.x

		a=pumpBeam.profile(esXPos)
		b=stedBeam.profile(esXPos)
		plt.plot(esXPos/1.0E-6, a)
		plt.plot(esXPos/1.0E-6, b)

		for idx in range(electronicSystems.N):
			plt.plot(esXPos[idx]/1.0E-6, 0.0009, 'ro')		
			if electronicSystems.isPopulated(idx):
				plt.plot(esXPos[idx]/1.0E-6, 0.001, 'g^')

		plt.show()
		

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
