
from Utility import EvolutionRecorder
import numpy as np

#==============================================================================
class ElectronicSystem(object):
	#--------------------------------------------------------------------------
	def __init__(self, N, REidx, xPos, pumpBeam, stedBeam):
		"""creates N electronic systems. RE must be specified separately in here"""

		self._pumpBeam = pumpBeam
		self._stedBeam = stedBeam

		self._numberElements = N
		self._rareEarthIndex = REidx
		self._isPopulated = np.zeros((self._numberElements, self._numberElements), dtype=bool)
		self._isRareEarth = np.zeros((self._numberElements, self._numberElements), dtype=bool)

		self.x = np.array(xPos) # setting up x positions of all elements
		self.y = np.array(xPos) # setting y positions to be same like x.

		# init the rare earth in the electronic systems
		self._states = dict()
		self._states['ground']  = 1
		self._states['excited'] = 2
		self._states['ionized'] = 3

		self._isRareEarth[REidx] = True
		self._isPopulated[REidx] = True
		self._reState = self._states['ground']

		# some measures
		self._reGroundStateCounter  = 0
		self._reExcitedStateCounter = 0

		self.groundStateEvolution = EvolutionRecorder('ground state', 'sim step', 'N')
		self.excitedStateEvolution = EvolutionRecorder('excited state', 'sim step', 'N')

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	def setupTransitionProbabilities(self, gamma, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE):
		pumpIntensity = self._pumpBeam.profile(self.x)
		stedIntensity = self._stedBeam.profile(self.x)

		# for electron traps
		self._probIonizeET = pumpIntensity + stedIntensity
		self._probIonizeET /= np.max(self._probIonizeET)

		# for rare earth
		self._probExciteRE  = pumpIntensity * sigPumpRE
		self._probIonizeRE  = (pumpIntensity + stedIntensity) * sigIonizeRE
		self._probRepumpRE  = pumpIntensity * sigRepumpRE
		#self._probDecayRE   = np.array([gamma for i in range(self._numberElements)])
		self._probDecayRE   = np.full((self._numberElements, self._numberElements), gamma)
		self._probDepleteRE = stedIntensity * sigStedRE

		totProbRE = self._probExciteRE + self._probIonizeRE + self._probRepumpRE + self._probDecayRE + self._probDepleteRE

		maxTotProbRE = np.max(totProbRE)
		totProbRE /= maxTotProbRE
		self._probExciteRE  /= maxTotProbRE
		self._probIonizeRE  /= maxTotProbRE
		self._probRepumpRE  /= maxTotProbRE
		self._probDecayRE   /= maxTotProbRE
		self._probDepleteRE /= maxTotProbRE

		self._probDecayRE   = self._probDecayRE[self._rareEarthIndex]
		self._probIonizeRE  = self._probDecayRE + self._probIonizeRE[self._rareEarthIndex]
		self._probExciteRE  = self._probIonizeRE + self._probExciteRE[self._rareEarthIndex]
		self._probRepumpRE  = self._probExciteRE + self._probRepumpRE[self._rareEarthIndex]
		self._probDepleteRE = self._probRepumpRE + self._probDepleteRE[self._rareEarthIndex]

		totProbREval = 10.*(self._probExciteRE + self._probIonizeRE + self._probRepumpRE + self._probDecayRE + self._probDepleteRE)
		self._probExciteRE  /= totProbREval
		self._probIonizeRE  /= totProbREval
		self._probRepumpRE  /= totProbREval
		self._probDecayRE   /= totProbREval
		self._probDepleteRE /= totProbREval

		self._normProbability = self._probIonizeET
		self._normProbability[self._rareEarthIndex] = totProbRE[self._rareEarthIndex]

	#--------------------------------------------------------------------------
	@property
	def x(self):
		return self._xPos

	@x.setter
	def x(self, value):
		self._xPos = value

	#--------------------------------------------------------------------------
	@property
	def y(self):
		return self._yPos

	@x.setter
	def y(self, value):
		self._yPos = value

	#--------------------------------------------------------------------------
	@property
	def population(self):
		return self._isPopulated

	#--------------------------------------------------------------------------
	@property
	def N(self):
		return self._numberElements**2

	#--------------------------------------------------------------------------
	@property
	def reState(self):
		return self._reState

	@reState.setter
	def reState(self, value):
		self._reState = value

	#--------------------------------------------------------------------------
	@property
	def rareEarthGroundStateCounter(self):
		return self._reGroundStateCounter

	#--------------------------------------------------------------------------
	@property
	def rareEarthExcitedStateCounter(self):
		return self._reExcitedStateCounter

	#--------------------------------------------------------------------------
	def resetEvolutionCounters(self):
		self._reGroundStateCounter  = 0
		self._reExcitedStateCounter = 0

	#--------------------------------------------------------------------------
	def recordREstate(self):
		if self._reState == self._states['ground']:
			self._reGroundStateCounter += 1

		elif self._reState == self._states['excited']:
			self._reExcitedStateCounter += 1

		else:
			pass

	#--------------------------------------------------------------------------
	def getPosition(self, idx):
		#print idx
		return np.array([self._xPos[idx[0]], self._yPos[idx[1]]])

	#--------------------------------------------------------------------------
	def isPopulated(self, idx):
		return self._isPopulated[idx]

	#--------------------------------------------------------------------------
	def isRareEarth(self, idx):
		return self._isRareEarth[idx]

	#--------------------------------------------------------------------------
	def restore(self):
		if self.reState == self._states['ionized']:
			self._isPopulated[self._rareEarthIndex] = True
			self.reState = self._states['ground']
			#print "restored"
			return 2
		else:
			#print "cannot restore RE, not in ionized state"
			return 0		

	#--------------------------------------------------------------------------
	def excite(self):
		if self.reState == self._states['ground']:
			self.reState = self._states['excited']
			#print "excited"
		else:
			#print "cannot excite RE, not in ground state"
			pass
	
		return 0

	#--------------------------------------------------------------------------
	def ionize(self, idx):
		if self.isPopulated(idx):
			if self.isRareEarth(idx):
				if self.reState == self._states['excited']:
					self.reState = self._states['ionized']
					self._isPopulated[idx] = False
					#print "ionized RE"
					return 1
				else:
					#print "cannot ionize RE, not in excited state"
					return 0
			else:
				self._isPopulated[idx] = False
				#print "ionized ET"
				return 1

		else:
			return 0

	#--------------------------------------------------------------------------
	def decay(self):
		if self.reState == self._states['excited']:
			#print "decay"
			self.setElectronToGroundState()

		return 0

	#--------------------------------------------------------------------------
	def deplete(self):
		if self.reState == self._states['excited']:
			#print "depleted"
			self.setElectronToGroundState()

		return 0

	#--------------------------------------------------------------------------
	def recombine(self, idx):
		if not self.isPopulated(idx):
			self._isPopulated[idx] = True
			#print "recombined ET"

			if self.isRareEarth(idx):
				self.reState = self._states['excited']
				#print "recombined RE"

			return 0

	#--------------------------------------------------------------------------
	def setElectronToGroundState(self):
		self.reState = self._states['ground']

	#--------------------------------------------------------------------------
	def actOnElectronTrap(self, idx, probability):
		if probability <= self._normProbability[idx]:
			return self.ionize(idx)

		else:
			return 0

	#--------------------------------------------------------------------------
	def actOnRareEarth(self, probability):
		#print probability
		if probability <= self._probDecayRE:
			#print "try to decay"
			return self.decay()

		elif probability <= self._probIonizeRE:
			#print "try to ionize"
			return self.ionize(self._rareEarthIndex)

		elif probability <= self._probExciteRE:
			#print "try to excite"
			return self.excite()

		elif probability <= self._probRepumpRE:
			#print "try to repump"
			return self.restore()

		elif probability <= self._probDepleteRE:
			#print "try to deplete"
			return self.deplete()

		else:
			return 0
