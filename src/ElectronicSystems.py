
import copy

#==============================================================================
class Electron(object):
	#--------------------------------------------------------------------------
	def __init__(self, x=0.0, travelRange=25E-9):
		self._xPos = x
		self._travelRange = travelRange

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def x(self):
		return self._xPos

	@x.setter
	def x(self, value):
		self._xPos = value

	#--------------------------------------------------------------------------
	@property
	def travelRange(self):
		return self._travelRange

	@travelRange.setter
	def travelRange(self, value):
		self._travelRange = value

#==============================================================================
class ElectronicSystemBase(object):
	#--------------------------------------------------------------------------
	def __init__(self, x=0.0):
		self._xPos = x
		self._localLaserPower = 0.0
		self._localPumpPower  = 0.0
		self._localStedPower  = 0.0

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def isPopulated(self):
		return hasattr(self, '_electron')

	#--------------------------------------------------------------------------
	@property
	def x(self):
		return self._xPos

	@x.setter
	def x(self, value):
		self._xPos = value

	#--------------------------------------------------------------------------
	# represents the total local laser power (pump + sted)
	@property
	def localLaserPower(self):
		return self._localLaserPower

	#--------------------------------------------------------------------------
	# represents the total power of the sted beam
	@property
	def localStedPower(self):
		return self._localStedPower

	@localStedPower.setter
	def localStedPower(self, value):
		self._localStedPower = value
		self._localLaserPower = self.localPumpPower + value

	#--------------------------------------------------------------------------
	# represents the total power of the pump beam
	@property
	def localPumpPower(self):
		return self._localPumpPower

	@localPumpPower.setter
	def localPumpPower(self, value):
		self._localPumpPower = value
		self._localLaserPower = self.localStedPower + value

	#--------------------------------------------------------------------------
	# represents the probability of getting ionized by the STED beam
	@property
	def probabilityStedIonize(self):
		return self._probStedIonize

	@probabilityStedIonize.setter
	def probabilityStedIonize(self, value):
		self._probStedIonize = value

	#--------------------------------------------------------------------------
	# represents the probability of getting ionized by the pump beam
	@property
	def probabilityPumpIonize(self):
		return self._probPumpIonize

	@probabilityPumpIonize.setter
	def probabilityPumpIonize(self, value):
		self._probPumpIonize = value

	#--------------------------------------------------------------------------
	# represents the probability to recombine with an electron from the
	# conduction band
	@property
	def probabilityRecombine(self):
		return self._probRecombine

	@probabilityRecombine.setter
	def probabilityRecombine(self, value):
		self._probRecombine = value


#==============================================================================
class ElectronTrap(ElectronicSystemBase):
	#--------------------------------------------------------------------------
	def __init__(self, x=0.0):
		ElectronicSystemBase.__init__(self, x=x)

	#--------------------------------------------------------------------------
	def ionize(self):
		"""Check if electron is present. Returns the present electron
		   and deletes it afterwards locally (here in the electron trap)."""
		if self.isPopulated:
			tmp = self._electron # former deepcopy
			del self._electron
			return tmp
		else:
			return None

	#--------------------------------------------------------------------------
	def recombine(self, electron):
		"""Takes an electron from the conduction band to get populated."""
		self._electron = electron # former deepcopy
		self._electron.x = self.x

	#--------------------------------------------------------------------------
	def __str__(self):
		return "electron trap"


#==============================================================================
class RareEarth(ElectronicSystemBase):
	#--------------------------------------------------------------------------
	def __init__(self, x=0.0, electronTravelRange=25E-9):
		ElectronicSystemBase.__init__(self, x)

		self._states = dict()
		self._states['ground'] = 1
		self._states['excited'] = 2
		self._states['ionized'] = 3

		self._electron = Electron(x=self.x, travelRange=electronTravelRange)
		self.state = self._states['ground']

		#self.probabilityDecay = 0.01

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, value):
		self._state = value

	#--------------------------------------------------------------------------
	# represents the probability of the RE's excited state to decay to
	# it's ground state
	@property
	def probabilityDecay(self):
		return self._probDecay

	@probabilityDecay.setter
	def probabilityDecay(self, value):
		self._probDecay = value

	#--------------------------------------------------------------------------
	# represents the probability of the rare earth to get to excited state
	# by the pump beam
	@property
	def probabilityPumpExcite(self):
		return self._probPumpExcite

	@probabilityPumpExcite.setter
	def probabilityPumpExcite(self, value):
		self._probPumpExcite = value

	#--------------------------------------------------------------------------
	# represents the probability of the rare earth to get depleted by the
	# STED beam
	@property
	def probabilityStedDeplete(self):
		return self._probStedDeplete

	@probabilityStedDeplete.setter
	def probabilityStedDeplete(self, value):
		self._probStedDeplete = value

	#--------------------------------------------------------------------------
	# represents the probability of the rare earth to get restored with an
	# electron from the valence band to it's ground state with the pump beam.
	@property
	def probabilityPumpRestore(self):
		return self._probPumpRestore

	@probabilityPumpRestore.setter
	def probabilityPumpRestore(self, value):
		self._probPumpRestore = value

	#--------------------------------------------------------------------------
	def excite(self):
		"""Transition from ground state to excited state."""
		if self.isPopulated and self.state == self._states['ground']:
			self.state = self._states['excited']

	#--------------------------------------------------------------------------
	def ionize(self):
		if self.state == self._states['excited']:
			tmp = self._electron # former deepcopy
			del self._electron
			self.state = self._states['ionized']
			return tmp
		else:
			return None

	#--------------------------------------------------------------------------
	def recombine(self, electron):
		"""Recombine with an electron from the conduction band, which
		   will end up in excited state."""
		self._electron = electron # former deepcopy
		self._electron.x = self.x
		self.state = self._states['excited']

	#--------------------------------------------------------------------------
	def restore(self, electron):
		"""Recombine with an electron from the valence band, which
		   will end up in ground state."""
		if not self.isPopulated:
			self._electron = electron # former deepcopy
			self._electron.x = self.x
			self.state = self._states['ground']
		else:
			pass

	#--------------------------------------------------------------------------
	def decay(self):
		if self.state == self._states['excited']:
			self.state = self._states['ground']

	#--------------------------------------------------------------------------
	def deplete(self):
		self.decay()
		
	#--------------------------------------------------------------------------
	def __str__(self):
		return "rare earth"