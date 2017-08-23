
from Utility import EvolutionRecorder
import numpy as np

#==============================================================================
class ElectronicSystem(object):
	#--------------------------------------------------------------------------
	def __init__(self, RExPos, REyPos, ETxPos, ETyPos, pumpBeam, stedBeam):
		
		if ETxPos.size != ETyPos.size:
			raise ValueError("x and y position array for ET must be of same size.")

		if RExPos.size != REyPos.size:
			raise ValueError("x and y position array for RE must be of same size.")

		# access indices for an electronic systems (for both ET and RE)
		self.idx = dict()
		self.idx['x']           = 0
		self.idx['y']           = 1
		self.idx['z']           = 2
		self.idx['isRE']        = 3
		self.idx['isPopulated'] = 4
		self.idx['pDecay']      = 5
		self.idx['pIonize']     = 6
		self.idx['pExcite']     = 7
		self.idx['pRepump']     = 8
		self.idx['pDeplete']    = 9
		self.idx['reState']     = 10
		self.idx['groundStateCounter']  = 11
		self.idx['excitedStateCounter'] = 12

		# init the rare earth in the electronic systems
		self._states = dict()
		self._states['ground']  = 1.0
		self._states['excited'] = 2.0
		self._states['ionized'] = 3.0

		# setup and initialize electron traps
		self.electronTraps = np.zeros((ETxPos.size, len(self.idx.keys())))
		self.electronTraps[:, self.idx['x']] = np.array(ETxPos)
		self.electronTraps[:, self.idx['y']] = np.array(ETyPos)
		self.electronTraps[:, self.idx['z']] = np.zeros(ETxPos.size)

		# setup and initialize rare earths
		self.rareEarths = np.zeros((RExPos.size, len(self.idx.keys())))
		self.rareEarths[:, self.idx['x']] = np.array(RExPos)
		self.rareEarths[:, self.idx['y']] = np.array(REyPos)
		self.rareEarths[:, self.idx['z']] = np.zeros(self.rareEarths.shape[0])
		self.rareEarths[:, self.idx['isRE']] = np.ones(self.rareEarths.shape[0])
		self.rareEarths[:, self.idx['isPopulated']] = np.ones(self.rareEarths.shape[0])
		self.rareEarths[:, self.idx['reState']] = np.full(self.rareEarths.shape[0], self._states['ground'])

		# keep laser beams
		self._pumpBeam = pumpBeam
		self._stedBeam = stedBeam

	#--------------------------------------------------------------------------
	def setupTransitionProbabilities(self, gammaRE, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE):
		pumpIntensityET = self._pumpBeam.profile(self.electronTraps[:, self.idx['x']], self.electronTraps[:, self.idx['y']])
		stedIntensityET = self._stedBeam.profile(self.electronTraps[:, self.idx['x']], self.electronTraps[:, self.idx['y']])

		pumpIntensityRE = self._pumpBeam.profile(self.rareEarths[:, self.idx['x']], self.rareEarths[:, self.idx['y']])
		stedIntensityRE = self._stedBeam.profile(self.rareEarths[:, self.idx['x']], self.rareEarths[:, self.idx['y']])

		# probability for ionizing the electron traps
		pIonizeET = pumpIntensityET + stedIntensityET
		pIonizeET /= np.max(pIonizeET)

		self.electronTraps[:, self.idx['pIonize']] = pIonizeET

		# for rare earth
		probExciteRE  = pumpIntensityRE * sigPumpRE
		probIonizeRE  = (pumpIntensityRE + stedIntensityRE) * sigIonizeRE
		probRepumpRE  = pumpIntensityRE * sigRepumpRE
		probDecayRE   = np.full(self.rareEarths.shape[0], gammaRE)
		probDepleteRE = stedIntensityRE * sigStedRE

		totProbRE = probExciteRE + probIonizeRE + probRepumpRE + probDecayRE + probDepleteRE

		maxTotProbRE   = np.max(totProbRE)
		totProbRE     /= maxTotProbRE
		probExciteRE  /= maxTotProbRE
		probIonizeRE  /= maxTotProbRE
		probRepumpRE  /= maxTotProbRE
		probDecayRE   /= maxTotProbRE
		probDepleteRE /= maxTotProbRE

		#self._probDecayRE   = self._probDecayRE[self._rareEarthIndex]
		probIonizeRE  = probDecayRE  + probIonizeRE
		probExciteRE  = probIonizeRE + probExciteRE
		probRepumpRE  = probExciteRE + probRepumpRE
		probDepleteRE = probRepumpRE + probDepleteRE

		totProbREval   = 10.0 * (probExciteRE + probIonizeRE + probRepumpRE + probDecayRE + probDepleteRE)
		probExciteRE  /= totProbREval
		probIonizeRE  /= totProbREval
		probRepumpRE  /= totProbREval
		probDecayRE   /= totProbREval
		probDepleteRE /= totProbREval

		self.rareEarths[:, self.idx['pDecay']]   = probDecayRE
		self.rareEarths[:, self.idx['pIonize']]  = probIonizeRE
		self.rareEarths[:, self.idx['pExcite']]  = probExciteRE
		self.rareEarths[:, self.idx['pRepump']]  = probRepumpRE
		self.rareEarths[:, self.idx['pDeplete']] = probDepleteRE

		#self._normProbability = pIonizeET
		#self._normProbability[self._rareEarthIndex] = totProbRE[self._rareEarthIndex]

		# concatenate electron traps and rare earths to a single array system
		self.electronicSystem = np.vstack((self.electronTraps, self.rareEarths))

		self.resetRareEarthEvolutionCounters()

	#--------------------------------------------------------------------------
	@property
	def x(self):
		return self.electronicSystem[:, self.idx['x']]

	#--------------------------------------------------------------------------
	@property
	def y(self):
		return self.electronicSystem[:, self.idx['y']]

	#--------------------------------------------------------------------------
	@property
	def z(self):
		return self.electronicSystem[:, self.idx['z']]

	#--------------------------------------------------------------------------
	@property
	def N(self):
		"""Returns the total number of electronic systems.
		   That means electron traps and rare earths."""
		return self.electronicSystem.shape[0]

	#--------------------------------------------------------------------------
	@property
	def rareEarthIndices(self):
		"""Returns an array of integers, which represent the
		   indices of the rare earths in the electronic systems."""
		return np.where(self.electronicSystem[:, self.idx['isRE']] == 1.0)[0]

	#--------------------------------------------------------------------------
	@property
	def electronTrapIndices(self):
		"""Returns an array of integers, which represent the
		   indices of the electron traps in the electronic systems."""
		return np.where(self.electronicSystem[:, self.idx['isRE']] == 0.0)[0]

	#--------------------------------------------------------------------------
	@property
	def population(self):
		return self.electronicSystem[:, self.idx['isPopulated']]

	#--------------------------------------------------------------------------
	def rareEarthGroundStateCounter(self, idx):
		return self.electronicSystem[idx][self.idx['groundStateCounter']]

	#--------------------------------------------------------------------------
	def rareEarthExcitedStateCounter(self, idx):
		return self.electronicSystem[idx][self.idx['excitedStateCounter']]

	#--------------------------------------------------------------------------
	def resetRareEarthEvolutionCounters(self):
		for idx in self.rareEarthIndices:
			self.electronicSystem[idx][self.idx['groundStateCounter']] = 0.0
			self.electronicSystem[idx][self.idx['excitedStateCounter']] = 0.0

	#--------------------------------------------------------------------------
	def recordREstates(self):
		for idx in self.rareEarthIndices:
			if self.electronicSystem[idx][self.idx['reState']] == self._states['ground']:
				self.electronicSystem[idx][self.idx['groundStateCounter']] += 1.0

			elif self.electronicSystem[idx][self.idx['reState']] == self._states['excited']:
				self.electronicSystem[idx][self.idx['excitedStateCounter']] += 1.0

			else:
				pass

	#--------------------------------------------------------------------------
	def getPosition(self, idx):
		"""Returns the absolute position of an electronic system."""
		return self.electronicSystem[idx][:self.idx['z'] + 1]

	#--------------------------------------------------------------------------
	def isPopulated(self, idx):
		return self.electronicSystem[idx][self.idx['isPopulated']]

	#--------------------------------------------------------------------------
	def populate(self, idx):
		self.electronicSystem[idx][self.idx['isPopulated']] = 1.0

	#--------------------------------------------------------------------------
	def isRareEarth(self, idx):
		return self.electronicSystem[idx][self.idx['isRE']]

	#--------------------------------------------------------------------------
	@property
	def potentialRecombinationIndices(self):
		return np.where(self.electronicSystem[:, self.idx['isPopulated']] == 0.0)[0]

	#--------------------------------------------------------------------------
	def restore(self, idx):
		if self.electronicSystem[idx][self.idx['reState']] == self._states['ionized']:
			self.populate(idx)
			self.electronicSystem[idx][self.idx['reState']] = self._states['ground']
			return 2

		else:
			return 0

	#--------------------------------------------------------------------------
	def excite(self, idx):
		if self.electronicSystem[idx][self.idx['reState']] == self._states['ground']:
			self.electronicSystem[idx][self.idx['reState']] = self._states['excited']
	
		return 0

	#--------------------------------------------------------------------------
	def ionizeRE(self, idx):
		if self.electronicSystem[idx][self.idx['isPopulated']]:
			if self.electronicSystem[idx][self.idx['reState']] == self._states['excited']:
				self.electronicSystem[idx][self.idx['reState']] = self._states['ionized']
				self.electronicSystem[idx][self.idx['isPopulated']] = 0.0
				return 1

			else:
				return 0

		else:
			return 0

	#--------------------------------------------------------------------------
	def ionizeET(self, idx):
		if self.electronicSystem[idx][self.idx['isPopulated']]:
			self.electronicSystem[idx][self.idx['isPopulated']] = 0.0
			return 1

		else:
			return 0

	#--------------------------------------------------------------------------
	def decay(self, idx):
		if self.electronicSystem[idx][self.idx['reState']] == self._states['excited']:
			self.electronicSystem[idx][self.idx['reState']] = self._states['ground']

		return 0

	#--------------------------------------------------------------------------
	def deplete(self, idx):
		return self.decay(idx)

	#--------------------------------------------------------------------------
	def recombine(self, idx):
		if not self.isPopulated(idx):
			self.populate(idx)

			if self.isRareEarth(idx):
				self.electronicSystem[idx][self.idx['reState']] = self._states['excited']

			return 0

	#--------------------------------------------------------------------------
	def actOnElectronTrap(self, idx, probability):
		if probability <= self.electronicSystem[idx][self.idx['pIonize']]:
			return self.ionizeET(idx)

		else:
			return 0

	#--------------------------------------------------------------------------
	def actOnRareEarth(self, idx, probability):
		if probability <= self.electronicSystem[idx][self.idx['pDecay']]: #self._probDecayRE:
			return self.decay(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pIonize']]: #self._probIonizeRE:
			return self.ionizeRE(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pExcite']]: #self._probExciteRE:
			return self.excite(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pRepump']]: #self._probRepumpRE:
			return self.restore(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pDeplete']]: #self._probDepleteRE:
			return self.deplete(idx)

		else:
			return 0

	#--------------------------------------------------------------------------
	def actOnSystem(self, idx, probability):
		if self.electronicSystem[idx][self.idx['isRE']]:
			result = self.actOnRareEarth(idx, probability)
		else:
			result = self.actOnElectronTrap(idx, probability)

		return result