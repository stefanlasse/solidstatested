
import matplotlib.pyplot as plt

class EvolutionRecorder(object):
	#--------------------------------------------------------------------------
	def __init__(self, name='untitled', xlabel='x', ylabel='y'):
		self._x = list()
		self._y = list()
		self._name = name
		self._xlabel = xlabel
		self._ylabel = ylabel

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	def record(self, x, y):
		self._x.append(x)
		self._y.append(y)

	#--------------------------------------------------------------------------
	def plot(self):
		plt.plot(self._x, self._y)
		plt.title(self._name)
		plt.xlabel(self._xlabel)
		plt.ylabel(self._ylabel)
		plt.show()