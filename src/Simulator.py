
from ElectronicSystems import ElectronicSystem
from Crystal import ConductionBand, ValenceBand
from LaserProfiles import PumpBeam, StedBeam
from Visualizer import Visualizer
from Utility import EvolutionRecorder

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib as mpl
mpl.rcParams['font.size'] = 16
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import numpy as np
import sys

import threading

class SolidStateStedSimulator(threading.Thread):
	#--------------------------------------------------------------------------
	def __init__(self, nSimSteps=1E6, savePath="C:/"):
		"""
		Creates a simulator object for STED microscopy in solids.

		Parameters
		----------
		nSimSteps : int or float
			Number of iteration steps for the simulation
		savePath : string
			Path to where all simulation results are saved
		"""
		self.numberOfSimulationSteps = int(nSimSteps + 1)
		self.savePath = savePath

		threading.Thread.__init__(self)

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
		self.electronTrapXCoordinates = ETx
		self.electronTrapYCoordinates = ETy
		self.pumpAmplitude = pumpAmpl
		self.stedAmplitude = stedAmpl
		self.crossSections = cs
		self.electronTravelRange = eTR

		self.electronicSystemsPopulationDistribution = np.zeros(REx.size + ETx.size)

		# set up conduction and valence band
		self.vb = ValenceBand()

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
			if not simStep % progressUpdate:
				sys.stdout.write("\r%.0f %% "%(float(simStep)/float(self.numberOfSimulationSteps)*100.0))
				if int(float(simStep)/float(self.numberOfSimulationSteps)*100.0) == 99:
					print ""

			# build up list of random indices for electronic systems
			# to act on and always include all rare earths.
			self.randIndices[:numRandElectronicSystems] = np.random.choice(electronTrapIndices, numRandElectronicSystems, replace=False)
			self.randIndices[numRandElectronicSystems:] = rareEarthIndices
			np.random.shuffle(self.randIndices)

			for index in self.randIndices:
				randomNumber = np.random.rand()

				if index in rareEarthIndices: #if self.electronSystems.isRareEarth(index):
					result = self.electronSystems.actOnRareEarth(index, randomNumber)

					if result == 1:
						# RE got ionized, electron is absorbed by CB, now recombine to somewhere
						self.handleRecombination(self.electronSystems.getPosition(index))

					elif result == 2:
						# RE got repumped, catch electron from VB
						self.vb.donateElectron()
						#self.vb.evolutionDonatedElectrons.record(simStep, self.vb.donatedElectronCount)

				else:
					result = self.electronSystems.actOnElectronTrap(index, randomNumber)

					if result == 1:
						# ET got ionized, add electron to CB
						self.handleRecombination(self.electronSystems.getPosition(index))

			# record ground/excited state evolution (internal)
			self.electronSystems.recordREstates()

			# record ground/excited state evolution (binned)
			if not simStep % progressEvolutionRecord:
				for REidx in self.electronSystems.rareEarthIndices:
					self.evolutionRecoders[self.evolutionRecorderIdx[REidx]].record(simStep,
																					self.electronSystems.rareEarthGroundStateCounter(REidx),
																					self.electronSystems.rareEarthExcitedStateCounter(REidx))

				self.electronSystems.resetRareEarthEvolutionCounters()

			if not simStep % 2:
				self.recordElectronTrapPopulationDistribution(self.electronSystems.population)

		# after last simulation step
		#self.saveRareEarthPopulationEvolution()
		#self.saveElectronTrapPopulationDistribution()

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
	def saveRareEarthPopulationEvolution(self):
		plt.plot(self.evolutionRecoders[0]._t, self.evolutionRecoders[0]._g, label="ground")
		plt.plot(self.evolutionRecoders[0]._t, self.evolutionRecoders[0]._e, label="excited")
		plt.legend(loc='best')
		plt.title("RE pop evol, REpos=%s"%(repos.__str__()))
		plt.savefig('%spop_evol_REidx=%s_pump_%.2f_sted_%.0f.png'%(self.savePath, self.positionOfRareEarthCenter.__str__(), self.pumpAmplitude, self.stedAmplitude), dpi=300)
		plt.close()

	#--------------------------------------------------------------------------
	def saveElectronTrapPopulationDistribution(self):
		esXPos = self.electronSystems.x
		esYPos = self.electronSystems.y

		xv, yv = np.meshgrid(esXPos, esYPos)

		X = xv/1.0E-6
		Y = yv/1.0E-6
		Z = self.electronicSystemsPopulationDistribution
		Z[self.electronSystems._rareEarthIndex] = 0
		Z[Z < 1] = 1

		pumpProfile = self.pumpBeam.profile(esXPos, esYPos)
		stedProfile = self.stedBeam.profile(esXPos, esYPos)

		plt.pcolormesh(X, Y, Z, norm=LogNorm(vmin=Z.min(), vmax=Z.max()), cmap='viridis')
		plt.xlabel('x [um]')
		plt.ylabel('y [um]')
		repos = self.electronSystems.getPosition(self.electronSystems._rareEarthIndex)/1.0E-6
		plt.title("Trap dist, REpos=%s"%(repos.__str__()))

		plt.colorbar()
		plt.savefig('%sdistribution_REidx=%s_pump_%.2f_sted_%.0f.png'%(self.savePath, self.positionOfRareEarthCenter.__str__(), self.pumpAmplitude, self.stedAmplitude), dpi=300)
		plt.close()

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------

