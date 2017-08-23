

import numpy as np

#==============================================================================
class LaserBeamGaussian(object):
	#--------------------------------------------------------------------------
	def __init__(self, x=0.0, y=0.0, amplitude=1.0, wavelength=1E-6, numAperture=1.3):
		self._xPos = x
		self._yPos = y
		self._amplitude = amplitude
		self._wavelength = wavelength
		self._numAperture = numAperture

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
	def y(self):
		return self._yPos

	@y.setter
	def y(self, value):
		self._yPos = value

	#--------------------------------------------------------------------------
	@property
	def amplitude(self):
		return self._amplitude

	@amplitude.setter
	def amplitude(self, value):
		self._amplitude = value

	#--------------------------------------------------------------------------
	@property
	def wavelength(self):
		return self._wavelength

	@wavelength.setter
	def wavelength(self, value):
		self._wavelength = value

	#--------------------------------------------------------------------------
	@property
	def numericalAperture(self):
		return self._numAperture

	@numericalAperture.setter
	def numericalAperture(self, value):
		self._numAperture = value

	#--------------------------------------------------------------------------
	@property
	def fwhm(self):
		return self._getFWHM()

	#--------------------------------------------------------------------------
	@property
	def sigma(self):
		return self.fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))

	#--------------------------------------------------------------------------
	def getExponent(self, xVals, yVals):
		#yVals = yVals[:, np.newaxis]
		return 4.0*np.log(2.0) * ((np.square(xVals - self.x) + np.square(yVals - self.y))/np.square(self.fwhm))

	#--------------------------------------------------------------------------
	def _getFWHM(self):
		return self.wavelength/self.numericalAperture


#==============================================================================
class PumpBeam(LaserBeamGaussian):
	def profile(self, xVals, yVals):
		exponent = self.getExponent(xVals, yVals)
		return self.amplitude*np.exp(-exponent)

#==============================================================================
class StedBeam(LaserBeamGaussian):
	def profile(self, xVals, yVals):
		exponent = self.getExponent(xVals, yVals)
		return self.amplitude*exponent*np.exp(-exponent + 1.0)
