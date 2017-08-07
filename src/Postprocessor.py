

import glob
import os
import pickle
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['savefig.directory'] = os.chdir(os.path.dirname(__file__))
mpl.rcParams['font.size'] = 16

from lmfit.models import LorentzianModel, PowerLawModel, ConstantModel, ExpressionModel


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
			
			fitInit = model.eval(pars, x=x)
			out = model.fit(y, pars, x=x)

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


p = Postprocessor('../tmp/gamma_0.20_sigPumpRE_1.00_sigIonizeRE_5.00_sigRepumpRE_0.80_sigStedRE_0.50')
p.getAllResultFiles()
p.fit()