
from ElectronicSystems import ElectronicSystem
from Crystal import ConductionBand, ValenceBand
from LaserProfiles import PumpBeam, StedBeam
from Visualizer import Visualizer

import numpy as np
#import random
import copy
import sys

import threading

class SolidStateStedSimulator(threading.Thread):
	#--------------------------------------------------------------------------
	def __init__(self, nSimSteps=1E6, nET=300, eTR=25E-9, posRE=0, centerET=0.0, spanET=2E-6):
		self.numberOfSimulationSteps = int(nSimSteps)
		self.numberOfElectronTraps = nET
		self.electronTravelRange = eTR
		self.centerElectronTraps = centerET
		self.spanElectronTraps = spanET
		self.positionOfRareEarthCenter = int(posRE)

		self.visualizer = Visualizer()

		threading.Thread.__init__(self)

	#--------------------------------------------------------------------------
	def setupSimulation(self):

		np.random.seed()

		electronTrapXPosition = np.linspace(self.centerElectronTraps-self.spanElectronTraps/2.0,
											self.centerElectronTraps+self.spanElectronTraps/2.0,
											self.numberOfElectronTraps + 1)

		self.REindex = int((self.numberOfElectronTraps + 1)/2) + self.positionOfRareEarthCenter

		# set up conduction and valence band
		self.cb = ConductionBand(pos = electronTrapXPosition)
		self.vb = ValenceBand()

		# set up pump and sted beam
		self.pumpBeam = PumpBeam(x=0.0, amplitude=0.005, wavelength=470E-9, numAperture=1.3)
		self.stedBeam = StedBeam(x=0.0, amplitude=0.0,  wavelength=600E-9, numAperture=1.3)

		# set up electronic systems
		self.electronSystems = ElectronicSystem(N        = int(self.numberOfElectronTraps + 1),
			                                    REidx    = self.REindex,
			                                    xPos     = electronTrapXPosition,
			                                    pumpBeam = self.pumpBeam,
			                                    stedBeam = self.stedBeam)

		self.electronSystems.setupTransitionProbabilities(gamma       = 0.01, 
														  sigPumpRE   = 1.0,
														  sigIonizeRE = 1.0,
														  sigRepumpRE = 1.0,
														  sigStedRE   = 1.0)

	#--------------------------------------------------------------------------
	def run(self):

		rareEarthIndex = self.electronSystems._rareEarthIndex

		for simStep in xrange(self.numberOfSimulationSteps):
			if simStep % int(0.01*self.numberOfSimulationSteps) == 0:
				sys.stdout.write("\r%.0f %%"%(float(simStep)/float(self.numberOfSimulationSteps)*100.0))
				if int(float(simStep)/float(self.numberOfSimulationSteps)*100.0) == 99:
					print ""

			# record ground/excited state evolution
			if not simStep % 10000:
				self.electronSystems.groundStateEvolution.record(simStep, self.electronSystems.rareEarthGroundStateCounter)
				self.electronSystems.excitedStateEvolution.record(simStep, self.electronSystems.rareEarthExcitedStateCounter)
				# reset counters
				self.electronSystems._reGroundStateCounter = 0
				self.electronSystems._reExcitedStateCounter = 0

			randIndices = np.random.randint(low=0, high=self.numberOfElectronTraps+1, size=5)
			
			#TODO: this might slow whole system down
			if rareEarthIndex not in randIndices:
				randIndices = np.append(randIndices, rareEarthIndex)
				np.random.shuffle(randIndices)

			for index in randIndices:
				randomNumber = np.random.rand()

				if self.electronSystems.isRareEarth(index):
					result = self.electronSystems.actOnRareEarth(randomNumber)

					if result == 1:
						# RE got ionized, add electron to CB
						self.cb.absorbElectron(index)
						self.handleRecombination()

					elif result == 2:
						# RE got repumped, catch electron from VB
						self.vb.donateElectron()
						self.vb.evolutionDonatedElectrons.record(simStep, self.vb.donatedElectronCount)

				else:
					result = self.electronSystems.actOnElectronTrap(index, randomNumber)

					if result == 1:
						# ET got ionized, add electron to CB
						self.cb.absorbElectron(index)
						self.handleRecombination()

			# record how many electrons are in the CB every 1000 simulation steps
			#if simStep % int(0.01*self.numberOfSimulationSteps) == 0:
			#	self.cb.recordNumberOfAvailableElectrons(simStep)

	#--------------------------------------------------------------------------
	def handleRecombination(self):
		# now go through the conduction band's collected electrons and
		# find the ones which can recombine to either an electron trap
		# or to a rare earth.
		while self.cb.numberAvailableElectrons:
			self.electronIndices = self.cb.availableElectronIndices
			np.random.shuffle(self.electronIndices)
			electronIndex = self.electronIndices[0]
			self.possibleRecombinationSlots = np.where(self.electronSystems.population == False)[0]
			np.random.shuffle(self.possibleRecombinationSlots)
			for esIndex in self.possibleRecombinationSlots:
				esPos = self.electronSystems.getPosition(esIndex)
				electronPosition = self.cb.getElectronPosition(electronIndex)
				if np.abs(esPos - electronPosition) <= self.electronTravelRange:
						self.cb.donateElectron(electronIndex)
						self.electronSystems.recombine(esIndex)
						break

	#--------------------------------------------------------------------------
	def visualize(self):
		if np.sum(self.electronSystems.population) - self.vb.donatedElectronCount == 1:
			print "Number electrons is valid."
		else:
			print "Number electrons is INVALID."

		# plot electron distribution
		self.visualizer.visualize(self.pumpBeam, self.stedBeam, self.electronSystems, self.electronSystems._rareEarthIndex)

		# plot VB electron donation saturation
		self.vb.evolutionDonatedElectrons.plot()

		self.electronSystems.groundStateEvolution.plot()
		self.electronSystems.excitedStateEvolution.plot()

		



