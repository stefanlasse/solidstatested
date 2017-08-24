
import numpy as np

#==============================================================================
class ElectronicSystem(object):
	#--------------------------------------------------------------------------
	def __init__(self, RExPos, REyPos, ETxPos, ETyPos, pumpBeam, stedBeam):
		"""
		Construct a collection of electronic systems, which can consist of
		rare earth ions and electron traps.

		Parameters
		----------
		RExPos : array-like
			Represents the x-coordinates of the rare earths
		REyPos : array-like
			Represents the y-coordinates of the rare earths
		ETxPos : array-like
			Represents the x-coordinates of the electron traps
		ETyPos : array-like
			Represents the y-coordinates of the electron traps
		pumpBeam : Object of type LaserProfiles.PumpBeam()
			Represents the impact of the Gaussian-shaped excitation laser to the electronic systems
		stedBeam : Object of type LaserProfiles.StedBeam()
			Represents the impave of the donut-shaped STED laser
		"""
		
		if ETxPos.size != ETyPos.size:
			raise ValueError("x and y position array for ET must be of same size.")

		if RExPos.size != REyPos.size:
			raise ValueError("x and y position array for RE must be of same size.")

		# access indices for an electronic system (for both ET and RE)
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

		# encodings for the electronic states in a rare earth
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
		"""
		Builds and finalizes the collection of electronic systems and calculates
		the transition probabilities of each system depending on its position and
		assigned cross-section.

		Parameters
		----------
		gammaRE : float
			Spontaneous decay of the rare earth (no laser)
		sigPumpRE : float
			Cross-section of the rare earth for excitation (pump laser)
		sigIonizeRE : float
			Cross-section of the rare earth for ionization (pump and STED laser)
		sigRepumpRE : float
			Cross-section of the rare earth for restoration with electron from valence band (pump laser)
		sigStedRE : float
			Cross-section of the rare earth for depletion (STED laser)
		"""

		# calculate pump and STED laser intensities at each electron trap position
		pumpIntensityET = self._pumpBeam.profile(self.electronTraps[:, self.idx['x']], self.electronTraps[:, self.idx['y']])
		stedIntensityET = self._stedBeam.profile(self.electronTraps[:, self.idx['x']], self.electronTraps[:, self.idx['y']])

		# calculate pump and STED laser intensities at each rare earth position
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
		"""Returns an array of floats, which represents the x coordinates of
		all present electronic systems."""
		return self.electronicSystem[:, self.idx['x']]

	#--------------------------------------------------------------------------
	@property
	def y(self):
		"""Returns an array of floats, which represents the y coordinates of
		all present electronic systems."""
		return self.electronicSystem[:, self.idx['y']]

	#--------------------------------------------------------------------------
	@property
	def z(self):
		"""Returns an array of floats, which represents the z coordinates of
		all present electronic systems."""
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
		"""Retruns an array of floats, which represents the current
		population of all present electronic systems. A populated system
		is represented by 1.0, a non-populated system is represented by 0.0."""
		return self.electronicSystem[:, self.idx['isPopulated']]

	#--------------------------------------------------------------------------
	def rareEarthGroundStateCounter(self, idx):
		"""Returns a single float, which represents how often the electronic system
		with index idx (which must be a rare earth) has been in ground state since
		the last counter reset."""
		return self.electronicSystem[idx][self.idx['groundStateCounter']]

	#--------------------------------------------------------------------------
	def rareEarthExcitedStateCounter(self, idx):
		"""Returns a single float, which represents how often the electronic system
		with index idx (which must be a rare earth) has been in excited state since
		the last counter reset."""
		return self.electronicSystem[idx][self.idx['excitedStateCounter']]

	#--------------------------------------------------------------------------
	def resetRareEarthEvolutionCounters(self):
		"""Resets both the ground and excited state counter for all rare earths
		in the electronic system."""
		for idx in self.rareEarthIndices:
			self.electronicSystem[idx][self.idx['groundStateCounter']] = 0.0
			self.electronicSystem[idx][self.idx['excitedStateCounter']] = 0.0

	#--------------------------------------------------------------------------
	def recordREstates(self):
		"""Updates either the ground or excited state counter for all rare earths
		present in the system depending on their current electronic states."""
		for idx in self.rareEarthIndices:
			if self.electronicSystem[idx][self.idx['reState']] == self._states['ground']:
				self.electronicSystem[idx][self.idx['groundStateCounter']] += 1.0

			elif self.electronicSystem[idx][self.idx['reState']] == self._states['excited']:
				self.electronicSystem[idx][self.idx['excitedStateCounter']] += 1.0

			else:
				pass

	#--------------------------------------------------------------------------
	def getPosition(self, idx):
		"""Returns an array of floats, which represents the absolute position
		of the electronic system with index idx. The position is read as [x, y, z]."""
		return self.electronicSystem[idx][:self.idx['z'] + 1]

	#--------------------------------------------------------------------------
	def isPopulated(self, idx):
		"""Returns a single float, which indicated whether the electronic system
		with index idx is populated (1.0) or not populated (0.0)."""
		return self.electronicSystem[idx][self.idx['isPopulated']]

	#--------------------------------------------------------------------------
	def populate(self, idx):
		"""Populates the electronic system with index idx. This function
		doesn't take care whether this is physically possible or not, nor
		does it take care if the system to be populated is a rare earth or
		and must be put into some electronic state."""
		self.electronicSystem[idx][self.idx['isPopulated']] = 1.0

	#--------------------------------------------------------------------------
	def isRareEarth(self, idx):
		"""Returns a float representing whether the electronic system
		with index idx is a rare earth (1.0) or an electron trap (0.0)."""
		return self.electronicSystem[idx][self.idx['isRE']]

	#--------------------------------------------------------------------------
	@property
	def potentialRecombinationIndices(self):
		"""Retruns an array of integers representing the indices of all
		electronic systems which are not populated."""
		return np.where(self.electronicSystem[:, self.idx['isPopulated']] == 0.0)[0]

	#--------------------------------------------------------------------------
	def restore(self, idx):
		"""Populates a rare earth with index idx if it is currently ionized
		and puts its electron to ground state. Here, the electron comes from
		the valence band."""
		if self.electronicSystem[idx][self.idx['reState']] == self._states['ionized']:
			self.populate(idx)
			self.electronicSystem[idx][self.idx['reState']] = self._states['ground']
			return 2

		else:
			return 0

	#--------------------------------------------------------------------------
	def excite(self, idx):
		"""Puts the electron of the rare earth with index idx to excited state
		if it currently is in ground state."""
		if self.electronicSystem[idx][self.idx['reState']] == self._states['ground']:
			self.electronicSystem[idx][self.idx['reState']] = self._states['excited']
	
		return 0

	#--------------------------------------------------------------------------
	def ionizeRE(self, idx):
		"""Ionizes the rare earth with index idx if it currently is in excited state."""
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
		"""Ionizes an electron trap with index idx."""
		if self.electronicSystem[idx][self.idx['isPopulated']]:
			self.electronicSystem[idx][self.idx['isPopulated']] = 0.0
			return 1

		else:
			return 0

	#--------------------------------------------------------------------------
	def decay(self, idx):
		"""Puts the electron of a rare earth with index idx to ground state
		if it currently is in excited state."""
		if self.electronicSystem[idx][self.idx['reState']] == self._states['excited']:
			self.electronicSystem[idx][self.idx['reState']] = self._states['ground']

		return 0

	#--------------------------------------------------------------------------
	def deplete(self, idx):
		"""Puts the electron of a rare earth with index idx to ground state
		if it currently is in excited state."""
		return self.decay(idx)

	#--------------------------------------------------------------------------
	def recombine(self, idx):
		"""Populates an electronic system with index idx if it currently is not
		populated. If this system is a rare earth, its state is set to excited
		state. Here, the electron comes from the condiction band."""
		if not self.isPopulated(idx):
			self.populate(idx)

			if self.isRareEarth(idx):
				self.electronicSystem[idx][self.idx['reState']] = self._states['excited']

			return 0

	#--------------------------------------------------------------------------
	def actOnElectronTrap(self, idx, probability):
		"""Decides whether an operation is perfomed on an electron trap
		depending on its transition probability."""
		if probability <= self.electronicSystem[idx][self.idx['pIonize']]:
			return self.ionizeET(idx)

		else:
			return 0

	#--------------------------------------------------------------------------
	def actOnRareEarth(self, idx, probability):
		"""Decides whether and which operation is performed on a rare earth
		depending on the corresponding transition probabilities."""
		if probability <= self.electronicSystem[idx][self.idx['pDecay']]:
			return self.decay(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pIonize']]:
			return self.ionizeRE(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pExcite']]:
			return self.excite(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pRepump']]:
			return self.restore(idx)

		elif probability <= self.electronicSystem[idx][self.idx['pDeplete']]:
			return self.deplete(idx)

		else:
			return 0
