
from Utility import EvolutionRecorder

import numpy as np


#==============================================================================
class ConductionBand(object):
	#--------------------------------------------------------------------------
	def __init__(self, pos):
		self._electronPositionsX = pos[0]
		self._electronPositionsY = pos[1]
		self._isPopulated = np.zeros((self._electronPositionsX.size, self._electronPositionsY.size), dtype=bool)

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def availableElectronIndices(self):
		"""Returns array of indices for populated electron positions."""
		return np.vstack(np.where(self._isPopulated == True)).T

	#--------------------------------------------------------------------------
	@property
	def numberAvailableElectrons(self):
		"""Returns the amount of electrons in the conduction band."""
		return np.sum(self._isPopulated)

	#--------------------------------------------------------------------------
	def getElectronPosition(self, idx):
		"""Returns the position of a certain electron in the conduction band."""
		return np.array([self._electronPositionsX[idx[0]], self._electronPositionsY[idx[1]]])

	#--------------------------------------------------------------------------
	def absorbElectron(self, idx):
		"""Absorbs an electron from an electron trap or a rare earth."""
		self._isPopulated[idx] = True
		
	#--------------------------------------------------------------------------
	def donateElectron(self, idx):
		"""Releases a certain electron from the CB to an electron trap
		   or to a rare earth."""
		self._isPopulated[idx] = False


#==============================================================================
class ValenceBand(object):
	#--------------------------------------------------------------------------
	def __init__(self):
		self._donatedElectronCount = 0
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

