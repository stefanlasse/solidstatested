

import glob
import os
import pickle
import numpy as np

from lmfit.models import LorentzianModel, PowerLawModel
import matplotlib.pyplot as plt


class Postprocessor(object):
	#--------------------------------------------------------------------------
	def __init__(self):
		pass

	#--------------------------------------------------------------------------
	def __del__(self):
		pass

	#--------------------------------------------------------------------------
	def getAllResultFiles(self, directory='../tmp'):
		files = glob.glob(directory + '/' + '*.pys')
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

			#plt.plot(x, y, 'b-', label='data')
			#plt.plot(x, fitInit, 'k--', label='guess')
			#plt.plot(x, out.best_fit, 'r-', label='fit')
			#plt.title("pump = %.2f, sted = %.2f"%(pumpAmplitude, stedAmplitude))
			#plt.legend(loc='best')
			#plt.show()

			self.stedPowers.append(stedAmplitude)
			self.fwhm.append(out.params['l1_fwhm'].value)
		
		self.stedPowers, self.fwhm = zip(*sorted(zip(self.stedPowers, self.fwhm)))
		x = np.array(self.stedPowers)
		y = np.array(self.fwhm)
		x = x[np.where((x >= 5.0) & (x <= 40.0))]
		y = y[np.where((x >= 5.0) & (x <= 40.0))]

		model = PowerLawModel(prefix='pl1_')
		pars = model.guess(y, x=x)
		pars['pl1_exponent'].set(-2.0, min=-5.0, max=0.0)
		modEval = model.eval(pars, x=x)
		out = model.fit(y, pars, x=x)

		plt.plot(x, y, label='data')
		plt.plot(x, modEval, label='guess')
		plt.plot(x, out.best_fit, label='fit')
		plt.xlabel('STED intensity')
		plt.ylabel('FWHM')
		plt.title('exponent = %f'%out.params['pl1_exponent'])
		plt.legend(loc='best')
		plt.show()



	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------


p = Postprocessor()
p.getAllResultFiles()
p.fit()