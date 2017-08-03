
from Utility import EvolutionRecorder

import numpy as np
import matplotlib.pyplot as plt


#==============================================================================
class ConductionBand(object):
	#--------------------------------------------------------------------------
	def __init__(self, pos):
		self._electronPositionsX = pos[0]
		self._electronPositionsY = pos[1]
		self._isPopulated = np.zeros((self._electronPositionsX.size, self._electronPositionsY.size), dtype=bool)
		self._availableElectrons = 0

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def availableElectronIndices(self):
		"""Returns array of indices for populated electron positions"""
		return np.vstack(np.where(self._isPopulated == False)).T

	#--------------------------------------------------------------------------
	@property
	def numberAvailableElectrons(self):
		return np.sum(self._isPopulated)

	#--------------------------------------------------------------------------
	def getElectronPosition(self, idx):
		"""Returns the position of a certain electron in the conduction band."""
		return np.array([self._electronPositionsX[idx[0]], self._electronPositionsX[idx[0]]])

	#--------------------------------------------------------------------------
	def absorbElectron(self, idx):
		"""Absorbs an electron from an electron trap or a rare earth."""
		self._isPopulated[idx] = True
		
	#--------------------------------------------------------------------------
	def donateElectron(self, idx):
		"""Releases a certain electron from the CB to an electron trap
		   or to a rare earth."""
		self._isPopulated[idx] = False
		
	#--------------------------------------------------------------------------
	def recordNumberOfAvailableElectrons(self, simstep):
		pass



#==============================================================================
class ValenceBand(object):
	#--------------------------------------------------------------------------
	def __init__(self):
		self._donatedElectronCount = 0
		self._evolutionTime = list()
		self.evolutionDonatedElectrons = EvolutionRecorder('donated VB electrons', 'sim step', 'N')

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
	#def recordEvolution(self, simulationStep):
	#	self._evolutionTime.append(simulationStep)
	#	self._evolutionDonatedElectrons.append(self.donatedElectronCount)

	#--------------------------------------------------------------------------
	def plotElectronEvolution(self):
		plt.plot(self._evolutionTime, self._evolutionDonatedElectrons)
		plt.show()