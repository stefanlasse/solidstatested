
import matplotlib as mpl
mpl.rcParams['font.size'] = 16

import matplotlib.pyplot as plt

class EvolutionRecorder(object):
	#--------------------------------------------------------------------------
	def __init__(self, name='untitled', xlabel='x', ylabel='y'):
		self._t = list()
		self._g = list()
		self._e = list()
		self._name = name
		self._xlabel = xlabel
		self._ylabel = ylabel

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	def record(self, t, g, e):
		self._t.append(t)
		self._g.append(g)
		self._e.append(e)

	#--------------------------------------------------------------------------
	def plot(self):
		plt.plot(self._t, self._g, label="ground")
		plt.plot(self._t, self._e, label="excited")
		plt.title(self._name)
		plt.xlabel(self._xlabel)
		plt.ylabel(self._ylabel)
		plt.legend(loc='best')
		plt.show()
