
from ElectronicSystems import ElectronicSystem
#from Crystal import ConductionBand, ValenceBand
from LaserProfiles import PumpBeam, StedBeam
from Utility import EvolutionRecorder


import numpy as np
import sys

from multiprocessing import Process

class SolidStateStedSimulator(Process):
	#--------------------------------------------------------------------------
	def __init__(self, nSimSteps, resultContainer):
		"""
		Creates a simulator object for STED microscopy in solids.

		Parameters
		----------
		nSimSteps : int or float
			Number of iteration steps for the simulation
		resultContainer : multiprocessing.Queue
			Container to save simulation results
		"""
		super(SolidStateStedSimulator, self).__init__()

		self.numberOfSimulationSteps = int(nSimSteps + 1)
		self.resultContainer = resultContainer

	#--------------------------------------------------------------------------
	def setupSimulation(self, REx, REy, ETx, ETy, pumpAmpl=0.05, stedAmpl=0.5, cs=[1,1,1,1,1], eTR=25E-9):
		"""
		Configures all necessary parameters and creates the objects needed for the simulator.

		Parameters
		----------
		REx : array-like
			Represents the x-coordinates of the rare earths
		REy : array-like
			Represents the y-coordinates of the rare earths
		ETx : array-like
			Represents the x-coordinates of the electron traps
		ETy : array-like
			Represents the y-coordinates of the electron traps
		pumpAmpl : float
			Amplitude for the excitation laser beam
		stedAmpl : float
			Amplitude for the STED laser beam
		cs : array-like
			Cross-section for the rare earth in the order
			[gammaRE, sigPumpRE, sigIonizeRE, sigRepumpRE, sigStedRE]
		eTR : float
			Electron travel range
		"""
		self.rareEarthXCoordinates = np.array(REx)
		self.rareEarthYCoordinates = np.array(REy)
		self.electronTrapXCoordinates = np.array(ETx)
		self.electronTrapYCoordinates = np.array(ETy)
		self.pumpAmplitude = pumpAmpl
		self.stedAmplitude = stedAmpl
		self.crossSections = cs
		self.electronTravelRange = eTR

		self.electronicSystemsPopulationDistribution = np.zeros(REx.size + ETx.size)

		# set up valence band
		#self.vb = ValenceBand()

		# set up pump and sted beam
		self.pumpBeam = PumpBeam(x=0.0, y=0.0, amplitude=self.pumpAmplitude, wavelength=470E-9, numAperture=1.3)
		self.stedBeam = StedBeam(x=0.0, y=0.0, amplitude=self.stedAmplitude, wavelength=600E-9, numAperture=1.3)

		# set up electronic systems
		self.electronSystems = ElectronicSystem(RExPos   = self.rareEarthXCoordinates,
												REyPos   = self.rareEarthYCoordinates,
			                                    ETxPos   = self.electronTrapXCoordinates,
			                                    ETyPos   = self.electronTrapYCoordinates,
			                                    pumpBeam = self.pumpBeam,
			                                    stedBeam = self.stedBeam)

		self.electronSystems.setupTransitionProbabilities(gammaRE     = self.crossSections[0],
														  sigPumpRE   = self.crossSections[1],
														  sigIonizeRE = self.crossSections[2],
														  sigRepumpRE = self.crossSections[3],
														  sigStedRE   = self.crossSections[4])

		# set up evolution recorder for every rare earth in the system
		self.evolutionRecoders = list()
		self.evolutionRecorderIdx = dict()
		for reCnt, reIdx in zip(range(self.electronSystems.rareEarthIndices.size), self.electronSystems.rareEarthIndices):
			self.evolutionRecorderIdx[reIdx] = reCnt
			rePos = self.electronSystems.getPosition(reIdx)
			evRec = EvolutionRecorder('REpos=[%.2g, %.2g, %.2g]'%(rePos[0], rePos[1], rePos[2]), 'sim step', 'N')
			self.evolutionRecoders.append(evRec)

		np.random.seed()

	#--------------------------------------------------------------------------
	def run(self):
		rareEarthIndices = self.electronSystems.rareEarthIndices
		electronTrapIndices = self.electronSystems.electronTrapIndices

		# randomly choose ~1% of the available electron traps
		# to be handeled in a single simulation step
		numRandElectronicSystems = int(0.01 * electronTrapIndices.size + 1)

		# prepare array for random indices to loop over
		self.randIndices = np.zeros(numRandElectronicSystems + rareEarthIndices.size, dtype=np.int32)

		# some constants to check during simulation
		progressUpdate = int(0.01*self.numberOfSimulationSteps)
		progressEvolutionRecord = int(0.05*self.numberOfSimulationSteps)

		for simStep in xrange(self.numberOfSimulationSteps):
			#if not simStep % progressUpdate:
			#	sys.stdout.write("\r%.0f %% "%(float(simStep)/float(self.numberOfSimulationSteps)*100.0))
			#	if int(float(simStep)/float(self.numberOfSimulationSteps)*100.0) == 99:
			#		print ""

			# build up list of random indices for electronic systems
			# to act on and always include all rare earths.
			self.randIndices[:numRandElectronicSystems] = np.random.choice(electronTrapIndices, numRandElectronicSystems, replace=False)
			self.randIndices[numRandElectronicSystems:] = rareEarthIndices
			np.random.shuffle(self.randIndices)

			for index in self.randIndices:
				randomNumber = np.random.rand()

				if index in rareEarthIndices:
					result = self.electronSystems.actOnRareEarth(index, randomNumber)

					if result == 1:
						# RE got ionized, electron is absorbed by CB, now recombine to somewhere
						self.handleRecombination(self.electronSystems.getPosition(index))

					elif result == 2:
						# RE got repumped, catch electron from VB
						#self.vb.donateElectron()
						#self.vb.evolutionDonatedElectrons.record(simStep, self.vb.donatedElectronCount)
						pass

				else:
					result = self.electronSystems.actOnElectronTrap(index, randomNumber)

					if result == 1:
						# ET got ionized, add electron to CB and recombine to somewhere
						self.handleRecombination(self.electronSystems.getPosition(index))

			# record ground/excited state evolution (internal)
			self.electronSystems.recordREstates()

			# record ground/excited state evolution (binned for result)
			if not simStep % progressEvolutionRecord:
				for REidx in self.electronSystems.rareEarthIndices:
					self.evolutionRecoders[self.evolutionRecorderIdx[REidx]].record(simStep,
																					self.electronSystems.rareEarthGroundStateCounter(REidx),
																					self.electronSystems.rareEarthExcitedStateCounter(REidx))

				self.electronSystems.resetRareEarthEvolutionCounters()

			if not simStep % 2:
				self.recordElectronTrapPopulationDistribution(self.electronSystems.population)

		# after last simulation step
		self.finalize()

	#--------------------------------------------------------------------------
	def handleRecombination(self, electronPosition):
		# now go through the conduction band's collected electrons and
		# find the ones which can recombine to either an electron trap
		# or to a rare earth.
		self.possibleRecombinationSlots = self.electronSystems.potentialRecombinationIndices
		np.random.shuffle(self.possibleRecombinationSlots)

		for esIndex in self.possibleRecombinationSlots:
			esPosition = self.electronSystems.getPosition(esIndex)

			if np.linalg.norm(esPosition - electronPosition) <= self.electronTravelRange:
				self.electronSystems.recombine(esIndex)
				break

	#--------------------------------------------------------------------------
	def recordElectronTrapPopulationDistribution(self, distribution):
		self.electronicSystemsPopulationDistribution += distribution

	#--------------------------------------------------------------------------
	def finalize(self):
		# generate the values which are needed for further processing
		self.groundStateAverage = np.average(np.array_split(self.evolutionRecoders[0]._g, 2)[1])
		self.excitedStateAverage = np.average(np.array_split(self.evolutionRecoders[0]._e, 2)[1])

		# collect single results in a dictionary
		result = dict()
		result["reXpos"] = self.rareEarthXCoordinates
		result["reYpos"] = self.rareEarthYCoordinates
		result["groundStateAverage"] = self.groundStateAverage
		result["excitedStateAverage"] = self.excitedStateAverage
		result["rePopulationEvolution_time"] = self.evolutionRecoders[0]._t
		result["rePopulationEvolution_groundState"] = self.evolutionRecoders[0]._g
		result["rePopulationEvolution_excitedState"] = self.evolutionRecoders[0]._e
		result["pumpAmplitude"] = self.pumpAmplitude
		result["stedAmplitude"] = self.stedAmplitude
		result["crossSections"] = self.crossSections
		result["electronTravelRange"] = self.electronTravelRange
		result["electronTrapXCoordinates"] = self.electronTrapXCoordinates
		result["electronTrapYCoordinates"] = self.electronTrapYCoordinates
		result["populationDistribution"] = self.electronicSystemsPopulationDistribution

		self.resultContainer.put(result)
