

import glob
import os
import pickle
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['savefig.directory'] = os.chdir(os.path.dirname(__file__))
mpl.rcParams['font.size'] = 16

from lmfit.models import LorentzianModel, PowerLawModel, ConstantModel, ExpressionModel


#from mpl_toolkits.axes_grid1 import host_subplot
#import mpl_toolkits.axisartist as AA
#import matplotlib as mpl
#mpl.rcParams['font.size'] = 16
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm


	#--------------------------------------------------------------------------
	#def saveRareEarthPopulationEvolution(self):
	#	ax = plt.subplot(111)
	#	plt.plot(self.evolutionRecoders[0]._t, self.evolutionRecoders[0]._g, label="ground")
	#	plt.plot(self.evolutionRecoders[0]._t, self.evolutionRecoders[0]._e, label="excited")
	#	plt.xlabel("Time [sim steps]")
	#	plt.ylabel("N [1]")
	#	plt.legend(loc='best')
	#	plt.suptitle("RE evol, pos=[%.3g, %.3g]"%(self.rareEarthXCoordinates, self.rareEarthYCoordinates))
	#	plt.title("pump=%.2f , sted=%.2f"%(self.pumpAmplitude, self.stedAmplitude))
	#	return ax

	#--------------------------------------------------------------------------
	#def saveElectronTrapPopulationDistribution(self):
	#	esXPos = self.electronSystems.x
	#	esYPos = self.electronSystems.y

	#	xv, yv = np.meshgrid(esXPos, esYPos)

	#	X = xv/1.0E-6
	#	Y = yv/1.0E-6
	#	Z = self.electronicSystemsPopulationDistribution
	#	Z[self.electronSystems._rareEarthIndex] = 0
	#	Z[Z < 1.0] = 1.0

	#	plt.pcolormesh(X, Y, Z, norm=LogNorm(vmin=Z.min(), vmax=Z.max()), cmap='viridis')
	#	plt.xlabel('x [um]')
	#	plt.ylabel('y [um]')
	#	repos = self.electronSystems.getPosition(self.electronSystems._rareEarthIndex)/1.0E-6
	#	plt.title("Trap dist, REpos=%s"%(repos.__str__()))

	#	plt.colorbar()
	#	plt.savefig('%sdistribution_REidx=%s_pump_%.2f_sted_%.0f.png'%(self.savePath, self.positionOfRareEarthCenter.__str__(), self.pumpAmplitude, self.stedAmplitude), dpi=300)
	#	plt.close()















class Postprocessor(object):
	#--------------------------------------------------------------------------
	def __init__(self, directory):
		self.directory = directory

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	def getAllResultFiles(self):
		files = glob.glob(self.directory + '/' + '*.pys')
		self.resultfiles = list()
		for f in files:
			self.resultfiles.append(os.path.abspath(f))

	#--------------------------------------------------------------------------
	def fit(self):
		self.stedPowers = list()
		self.fwhm = list()

		for f in self.resultfiles:
			data = pickle.load(open(f, 'rb'))
			x = data['xPos']
			y = data['excitedStateAverage']
			pumpAmplitude = data['pumpAmplitude']
			stedAmplitude = data['stedAmplitude']

			if pumpAmplitude > 0.05 or pumpAmplitude < 0.05:
				continue

			model = LorentzianModel(prefix='l1_')
			pars = model.guess(y, x=x)

			pars['l1_center'].set(-1.5E-6, min=-1.5E-6-1E-10, max=-1.5E-6+1E-10)

			fitInit = model.eval(pars, x=x)
			out = model.fit(y, pars, x=x)

			print "%.3g, %.5g"%(out.params['l1_center'].value, out.params['l1_fwhm'].value)

			plt.plot(x, y, 'b-', label='data')
			plt.plot(x, fitInit, 'k--', label='guess')
			plt.plot(x, out.best_fit, 'r-', label='fit')
			plt.title("pump = %.2f, sted = %.2f"%(pumpAmplitude, stedAmplitude))
			plt.legend(loc='best')
			plt.show()

			self.stedPowers.append(stedAmplitude)
			self.fwhm.append(out.params['l1_fwhm'].value)

		self.stedPowers, self.fwhm = zip(*sorted(zip(self.stedPowers, self.fwhm)))
		x = np.array(self.stedPowers)
		y = np.array(self.fwhm)
		#x = x[np.where((x >= 5.0) & (x <= 40.0))]
		#y = y[np.where((x >= 5.0) & (x <= 40.0))]

		powerLaw = PowerLawModel(prefix='pl_')
		constant = ConstantModel(prefix='const_')
		pars = powerLaw.guess(y, x=x)
		pars.update(constant.make_params())
		pars['pl_exponent'].set(min=-10.0, max=0.0)
		pars['const_c'].set(0.0, min=0.0)
		model = powerLaw + constant

		#model = ExpressionModel('amp * 1.0/(x**exponent) + offset')
		#pars = model.make_params(amp=2.733E-7, exponent=0.5, offset=0.0)

		modEval = model.eval(pars, x=x)
		out = model.fit(y, pars, x=x)

		plt.rc('text', usetex=True)
		plt.plot(x, y, label='data')
		plt.plot(x, modEval, label='guess')
		plt.plot(x, out.best_fit, label='fit')
		plt.xlabel(r'$I_S$', fontsize=22)
		plt.ylabel(r'FWHM (PSF)', fontsize=22)
		plt.title(r'$\propto 1 / \sqrt[%.4f]{I_S}$'%(-1.0/out.params['pl_exponent'].value), fontsize=28)
		plt.legend(loc='best')
		plt.show()



	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------


p = Postprocessor('D:/STED_sim/2D/gamma_0.20_sigPumpRE_1.00_sigIonizeRE_15.00_sigRepumpRE_2.00_sigStedRE_1.00')
p.getAllResultFiles()
p.fit()
