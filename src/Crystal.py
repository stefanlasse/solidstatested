
import copy
import random
import matplotlib.pyplot as plt
from ElectronicSystems import Electron

#==============================================================================
class ConductionBand(object):
	#--------------------------------------------------------------------------
	def __init__(self):
		self._electrons = list()

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def availableElectrons(self):
		return len(self._electrons)

	#--------------------------------------------------------------------------
	def shuffleElectrons(self):
		random.shuffle(self._electrons)

	#--------------------------------------------------------------------------
	def getElectron(self, idx):
		return self._electrons[idx] # former deepcopy

	#--------------------------------------------------------------------------
	def absorbElectron(self, electron):
		"""Absorbs an electron from an electron trap or a rare earth."""
		if electron is not None:
			self._electrons.append(electron) # former deepcopy
		else:
			pass

	#--------------------------------------------------------------------------
	def donateElectron(self, idx):
		"""Releases a certain electron from the CB to an electron trap
		   or to a rare earth."""
		tmp = self._electrons[idx] # former deepcopy
		del self._electrons[idx]
		return tmp




#==============================================================================
class ValenceBand(object):
	#--------------------------------------------------------------------------
	def __init__(self):
		self._donatedElectronCount = 0
		self._evolutionTime = list()
		self._evolutionDonatedElectrons = list()

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def donatedElectronCount(self):
		return self._donatedElectronCount

	#--------------------------------------------------------------------------
	def donateElectron(self, x=0.0, travelRange=25E-9):
		self._donatedElectronCount += 1
		return Electron(x=x, travelRange=travelRange)

	#--------------------------------------------------------------------------
	def recordEvolution(self, simumationStep):
		self._evolutionTime.append(simumationStep)
		self._evolutionDonatedElectrons.append(self.donatedElectronCount)

	#--------------------------------------------------------------------------
	def plotElectronEvolution(self):
		plt.plot(self._evolutionTime, self._evolutionDonatedElectrons)
		plt.show()