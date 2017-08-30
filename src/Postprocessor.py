

import glob
import os
import pickle
import numpy as np


from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['savefig.directory'] = os.chdir(os.path.dirname(__file__))
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
def saveRareEarthPopulationEvolution(self):
	ax = plt.subplot(111)
	plt.plot(self.evolutionRecoders[0]._t, self.evolutionRecoders[0]._g, label="ground")
	plt.plot(self.evolutionRecoders[0]._t, self.evolutionRecoders[0]._e, label="excited")
	plt.xlabel("Time [sim steps]")
	plt.ylabel("N [1]")
	plt.legend(loc='best')
	plt.suptitle("RE evol, pos=[%.3g, %.3g]"%(self.rareEarthXCoordinates, self.rareEarthYCoordinates))
	plt.title("pump=%.2f , sted=%.2f"%(self.pumpAmplitude, self.stedAmplitude))
	return ax

#--------------------------------------------------------------------------
def saveElectronTrapPopulationDistribution(self):
	esXPos = self.electronSystems.x
	esYPos = self.electronSystems.y

	xv, yv = np.meshgrid(esXPos, esYPos)

	X = xv/1.0E-6
	Y = yv/1.0E-6
	Z = self.electronicSystemsPopulationDistribution
	Z[self.electronSystems._rareEarthIndex] = 0
	Z[Z < 1.0] = 1.0

	plt.pcolormesh(X, Y, Z, norm=LogNorm(vmin=Z.min(), vmax=Z.max()), cmap='viridis')
	plt.xlabel('x [um]')
	plt.ylabel('y [um]')
	repos = self.electronSystems.getPosition(self.electronSystems._rareEarthIndex)/1.0E-6
	plt.title("Trap dist, REpos=%s"%(repos.__str__()))

	plt.colorbar()
	plt.savefig('%sdistribution_REidx=%s_pump_%.2f_sted_%.0f.png'%(self.savePath, self.positionOfRareEarthCenter.__str__(), self.pumpAmplitude, self.stedAmplitude), dpi=300)
	plt.close()















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

		self.getData()

	#--------------------------------------------------------------------------
	def getData(self):
		self.data = dict()
		self.pumpAmplitudes = list()
		self.stedAmplitudes = list()

		for f in self.resultfiles:
			print f
			with open(f, 'rb') as fil:
				data = pickle.load(fil) # thats a list of dicts

			for i in range(len(data)):
				self.pumpAmplitudes.append(data[i]['pumpAmplitude'])
				self.stedAmplitudes.append(data[i]['stedAmplitude'])

		self.pumpAmplitudes = np.unique(self.pumpAmplitudes)
		self.stedAmplitudes = np.unique(self.stedAmplitudes)

		# build up data structure
		for pa in self.pumpAmplitudes:
			self.data[pa] = dict()
			for sa in self.stedAmplitudes:
				self.data[pa][sa] = [[],[]]

		# insert data to structure
		for f in self.resultfiles:
			with open(f, 'rb') as fil:
				data = pickle.load(fil) # thats a list of dicts

			for i in range(len(data)):
				pa = data[i]['pumpAmplitude']
				sa = data[i]['stedAmplitude']
				laserXpos = data[i]['laserXpos']
				esa = data[i]['excitedStateAverage']
				self.data[pa][sa][0].append(laserXpos)
				self.data[pa][sa][1].append(esa)

		# sort by laserXpos
		for pak in self.data.keys():
			for sak in self.data[pak].keys():
				laserXpos = np.array(self.data[pak][sak][0])
				esa = np.array(self.data[pak][sak][1])
				inds = laserXpos.argsort()
				self.data[pak][sak][0] = laserXpos[inds]
				self.data[pak][sak][1] = esa[inds]


	#--------------------------------------------------------------------------
	def fit(self, pa):

		self.stedPowers = list()
		self.fwhm = list()

		for sa in sorted(self.data[pa].keys()):	# ensure starting with lowest stedAmplitude
			x = self.data[pa][sa][0]
			y = self.data[pa][sa][1]/np.max(self.data[pa][sa][1])

			lorentzian = LorentzianModel(prefix='l1_')
			constant   = ConstantModel(prefix='const_')
			pars = lorentzian.guess(y, x=x)
			pars.update(constant.make_params())
			pars['const_c'].set(min=0.0)

			model = lorentzian + constant
			fitInit = model.eval(pars, x=x)
			out = model.fit(y, pars, x=x)

			#print "%.3g, %.5g"%(out.params['l1_center'].value, out.params['l1_fwhm'].value)

			#plt.plot(x, y, 'b-', label='data')
			#plt.plot(x, fitInit, 'k--', label='guess')
			#plt.plot(x, out.best_fit, 'r-', label='fit')
			#plt.title("pump=%.2f, sted=%.2f"%(pa, sa))
			#plt.legend(loc='best')
			#plt.show()

			self.stedPowers.append(sa)
			self.fwhm.append(out.params['l1_fwhm'].value)

		self.stedPowers, self.fwhm = zip(*sorted(zip(self.stedPowers, self.fwhm)))
		x = np.array(self.stedPowers)
		y = np.array(self.fwhm)

		powerLaw = PowerLawModel(prefix='pl_')
		constant = ConstantModel(prefix='const_')
		pars = powerLaw.guess(y, x=x)
		pars.update(constant.make_params())
		pars['pl_exponent'].set(min=-10.0, max=0.0)
		pars['const_c'].set(0.0, min=0.0)
		model = powerLaw + constant

		modEval = model.eval(pars, x=x)
		out = model.fit(y, pars, x=x)

		plt.plot(x, y, label='data')
		plt.plot(x, modEval, label='guess')
		plt.plot(x, out.best_fit, label='fit')
		plt.xlabel('I_STED', fontsize=22)
		plt.ylabel('FWHM (PSF)', fontsize=22)
		plt.title('pa=%.2f, exponent = %.4f'%(pa, out.params['pl_exponent'].value), fontsize=22)
		plt.legend(loc='best')
		plt.show()



	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------


p = Postprocessor('D:/STED_sim/new_2D/gamma_0.20_sigPumpRE_2.00_sigIonizeRE_10.00_sigRepumpRE_2.00_sigStedRE_0.00')
p.getAllResultFiles()

for pa in sorted(p.data.keys()):
	p.fit(pa)
