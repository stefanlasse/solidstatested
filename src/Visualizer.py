

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
	def visualize(self, pumpBeam, stedBeam, electronicSystems, reIndex, simStep):

		esXPos = electronicSystems.x

		a = pumpBeam.profile(esXPos)
		b = stedBeam.profile(esXPos)
		plt.plot(esXPos/1.0E-6, a/np.max(a))
		plt.plot(esXPos/1.0E-6, b/np.max(b))

		for idx in range(electronicSystems.N):
			if idx == reIndex:
				plt.plot(esXPos[idx]/1.0E-6, 0.1, 'bo')
			else:
				plt.plot(esXPos[idx]/1.0E-6, 0.1, 'ro')

			if electronicSystems.isPopulated(idx):
				plt.plot(esXPos[idx]/1.0E-6, 0.2, 'g^')


		plt.title("%07d"%simStep)

		path = "C:/Users/Stefan/Documents/projects/rate_equations/solidstatested/tmp/"
		plt.savefig("%scurr_state_REidx_%03d_pump_%.2f_sted_%.2f_cnt_%07d.png"%(path, reIndex, pumpBeam._amplitude, stedBeam._amplitude, simStep))
		plt.close()


		#plt.show()
		

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
