
from ElectronicSystems import ElectronTrap, RareEarth
from Crystal import ConductionBand, ValenceBand
from LaserProfiles import PumpBeam, StedBeam
from Visualizer import Visualizer

import numpy as np
import random
random.seed()
import copy
import sys

import threading

class SolidStateStedSimulator(threading.Thread):
	#--------------------------------------------------------------------------
	def __init__(self, nSimSteps=2E6, nET=300, eTR=25E-9, posRE=0.0, centerET=0.0, spanET=2E-6):
		self.numberOfSimulationSteps = int(nSimSteps)
		self.numberOfElectronTraps = nET
		self.electronTravelRange = eTR
		self.centerElectronTraps = centerET
		self.spanElectronTraps = spanET
		self.positionOfRareEarthCenter = posRE

		self.visualizer = Visualizer()

		threading.Thread.__init__(self)

	#--------------------------------------------------------------------------
	def setupSimulation(self):
		# set up conduction and valence band
		self.cb = ConductionBand()
		self.vb = ValenceBand()

		# set up pump and sted beam
		self.pumpBeam = PumpBeam(x=0.0, amplitude=0.005, wavelength=470E-9, numAperture=1.3)
		self.stedBeam = StedBeam(x=0.0, amplitude=0.01,  wavelength=600E-9, numAperture=1.3)

		electronTrapXPosition = np.linspace(self.centerElectronTraps-self.spanElectronTraps/2.0,
											self.centerElectronTraps+self.spanElectronTraps/2.0,
											self.numberOfElectronTraps)

		# set up electron traps
		self.electronTraps = [ElectronTrap(x=electronTrapXPosition[i]) for i in range(self.numberOfElectronTraps)]

		# set up rare earth
		self.rareEarth = RareEarth(x=self.positionOfRareEarthCenter)

		# calculate laser intensities at each electron trap and rare earth position
		for et in self.electronTraps:
			et.localPumpPower = self.pumpBeam.profile(et.x)
			et.localStedPower = self.stedBeam.profile(et.x)

		self.rareEarth.localPumpPower = self.pumpBeam.profile(self.rareEarth.x)
		self.rareEarth.localStedPower = self.stedBeam.profile(self.rareEarth.x)

		# set up transition probabilities for the electron traps
		for et in self.electronTraps:
			et.probabilityPumpIonize = 0.5*self.pumpBeam.profile(et.x)
			et.probabilityStedIonize = et.probabilityPumpIonize + 0.5*self.stedBeam.profile(et.x)
			#et.probabilityRecombine  = 0.0

		# set up transition probabilities for the rare earth
		self.rareEarth.probabilityDecay       = 0.01
		self.rareEarth.probabilityPumpIonize  = self.rareEarth.probabilityDecay       + self.pumpBeam.profile(self.rareEarth.x)
		self.rareEarth.probabilityStedIonize  = self.rareEarth.probabilityPumpIonize  + self.stedBeam.profile(self.rareEarth.x)
		self.rareEarth.probabilityPumpExcite  = self.rareEarth.probabilityStedIonize  + 0.2*self.pumpBeam.profile(self.rareEarth.x)
		self.rareEarth.probabilityPumpRestore = self.rareEarth.probabilityPumpExcite  + self.pumpBeam.profile(self.rareEarth.x)
		self.rareEarth.probabilityStedDeplete = self.rareEarth.probabilityPumpRestore + 0.2*self.stedBeam.profile(self.rareEarth.x)
		#rareEarth.probabilityRecombine   = 0.0

		self.electronTraps.append(copy.deepcopy(self.rareEarth))
		del self.rareEarth

	#--------------------------------------------------------------------------
	def run(self):
		for simStep in xrange(self.numberOfSimulationSteps):
			if simStep % int(0.05*self.numberOfSimulationSteps) == 0:
				sys.stdout.write("\r%.0f %%"%(float(simStep)/float(self.numberOfSimulationSteps)*100.0))
				random.seed()
				if int(float(simStep)/float(self.numberOfSimulationSteps)*100.0) == 95:
					print ""

			random.shuffle(self.electronTraps)

			for element in self.electronTraps:
				randomNumber = random.random()

				# handle rare earth
				if type(element) is RareEarth:
					if randomNumber <= element.probabilityDecay:
						element.decay()

					#elif randomNumber <= element.probabilityPumpIonize:
					#	self.cb.absorbElectron(element.ionize())

					elif randomNumber <= element.probabilityStedIonize:
						self.cb.absorbElectron(element.ionize())

					elif randomNumber <= element.probabilityPumpExcite:
						element.excite()

					elif randomNumber <= element.probabilityPumpRestore:
						if not element.isPopulated:
							element.restore(self.vb.donateElectron(element.x, self.electronTravelRange))
							self.vb.recordEvolution(simStep)

					elif randomNumber <= element.probabilityStedDeplete:
						element.deplete()

					else:
						pass

				# handle electron traps
				elif type(element) is ElectronTrap:
					if randomNumber <= element.probabilityStedIonize:
						self.cb.absorbElectron(element.ionize())

					else:
						pass

				else:
					pass

			#==========================================================================
			# recombination part
			#==========================================================================

			# record how many electrons are in the CB every 1000 simulation steps
			if simStep % int(0.01*self.numberOfSimulationSteps) == 0:
				self.cb.recordNumberOfAvailableElectrons(simStep)

			# now go through the conduction band's collected electrons and
			# find the ones which can recombine to either an electron trap
			# or to a rare earth.
			
			self.cb.shuffleElectrons()
			while self.cb.availableElectrons != 0:
				# get electron from shuffled list
				e = self.cb.getElectron(0)

				# search for an electron trap/RE which is within the electrons
				# travel range and where it can recombine if the trap/RE is currently ionized
				for element in self.electronTraps:
					if np.abs(element.x - e.x) <= e.travelRange and not element.isPopulated:
						element.recombine(self.cb.donateElectron(0))
						break

		#==========================================================================
		# here the simulation is finished, clean up now
		#==========================================================================

		# sort electronic systems by x-position
		self.electronTraps.sort(key=lambda element: element.x)

	#--------------------------------------------------------------------------
	def visualize(self):
		# plot electron distribution
		self.visualizer.visualize(self.pumpBeam, self.stedBeam, self.electronTraps)

		# plot VB electron donation saturation
		self.vb.plotElectronEvolution()

