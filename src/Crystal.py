

import numpy as np
import matplotlib.pyplot as plt


#==============================================================================
class ConductionBand(object):
	#--------------------------------------------------------------------------
	def __init__(self, pos):
		self._electronPositions = pos
		self._isPopulated = np.zeros(self._electronPositions.size, dtype=bool)
		self._availableElectrons = 0

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def availableElectrons(self):
		"""Returns array of indices for populated electron positions"""
		return np.nonzero(self._isPopulated)

	#--------------------------------------------------------------------------
	def getElectronPosition(self, idx):
		"""Returns the position of a certain electron in the conduction band."""
		return self._electronPositions[idx]

	#--------------------------------------------------------------------------
	def absorbElectron(self, idx):
		"""Absorbs an electron from an electron trap or a rare earth."""
		self._isPopulated[idx] = True
		self._availableElectrons += 1

	#--------------------------------------------------------------------------
	def donateElectron(self, idx):
		"""Releases a certain electron from the CB to an electron trap
		   or to a rare earth."""
		self._isPopulated[idx] = False
		self._availableElectrons -= 1

	#--------------------------------------------------------------------------
	def recordNumberOfAvailableElectrons(self, simstep):
		pass



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
	def donateElectron(self):
		self._donatedElectronCount += 1

	#--------------------------------------------------------------------------
	def recordEvolution(self, simulationStep):
		self._evolutionTime.append(simulationStep)
		self._evolutionDonatedElectrons.append(self.donatedElectronCount)

	#--------------------------------------------------------------------------
	def plotElectronEvolution(self):
		plt.plot(self._evolutionTime, self._evolutionDonatedElectrons)
		plt.show()