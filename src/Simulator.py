
from ElectronicSystems import ElectronicSystem
from Crystal import ConductionBand, ValenceBand
from LaserProfiles import PumpBeam, StedBeam
from Visualizer import Visualizer

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

import numpy as np
import sys

import threading

class SolidStateStedSimulator(threading.Thread):
	#--------------------------------------------------------------------------
	def __init__(self, nSimSteps=1E6, nET=300, eTR=25E-9, posRE=0, centerET=0.0, spanET=2E-6,
		               pumpAmpl=0.1, stedAmpl=0.1, satAvg=10000, savePath="C:/", cs=[0,0,0,0,0]):

		self.numberOfSimulationSteps = int(nSimSteps)
		self.numberOfElectronTraps = nET
		self.electronTravelRange = eTR
		self.centerElectronTraps = centerET
		self.spanElectronTraps = spanET
		self.positionOfRareEarthCenter = posRE
		self.pumpAmplitude = pumpAmpl
		self.stedAmplitude = stedAmpl
		self.saturationAveraging = int(satAvg)
		self.crossSections = cs

		self.savePath = savePath

		self.electronTrapPopulationDistribution = np.zeros((self.numberOfElectronTraps + 1,self.numberOfElectronTraps + 1))

		self.visualizer = Visualizer()
		threading.Thread.__init__(self)

	#--------------------------------------------------------------------------
	def setupSimulation(self):
		np.random.seed()

		electronTrapXPosition = np.linspace(self.centerElectronTraps-self.spanElectronTraps/2.0,
											self.centerElectronTraps+self.spanElectronTraps/2.0,
											self.numberOfElectronTraps + 1)
		electronTrapYPosition = np.copy(electronTrapXPosition)

		self.REindex = (int((self.numberOfElectronTraps + 1)/2) + self.positionOfRareEarthCenter,
						int((self.numberOfElectronTraps + 1)/2) + self.positionOfRareEarthCenter)
		# set up conduction and valence band
		self.cb = ConductionBand(pos=[electronTrapXPosition, electronTrapYPosition])
		self.vb = ValenceBand()

		# set up pump and sted beam
		self.pumpBeam = PumpBeam(x=0.0, amplitude=self.pumpAmplitude, wavelength=470E-9, numAperture=1.3)
		self.stedBeam = StedBeam(x=0.0, amplitude=self.stedAmplitude, wavelength=600E-9, numAperture=1.3)

		# set up electronic systems
		self.electronSystems = ElectronicSystem(N        = int(self.numberOfElectronTraps + 1),
			                                    REidx    = self.REindex,
			                                    xPos     = electronTrapXPosition,
			                                    pumpBeam = self.pumpBeam,
			                                    stedBeam = self.stedBeam)

		self.electronSystems.setupTransitionProbabilities(gamma       = self.crossSections[0], 
														  sigPumpRE   = self.crossSections[1],
														  sigIonizeRE = self.crossSections[2],
														  sigRepumpRE = self.crossSections[3],
														  sigStedRE   = self.crossSections[4])

		self.electronSystems.resetEvolutionCounters()

	#--------------------------------------------------------------------------
	def run(self):
		rareEarthIndex = self.electronSystems._rareEarthIndex

		for simStep in xrange(self.numberOfSimulationSteps):
			if simStep % int(0.1*self.numberOfSimulationSteps) == 0:
				sys.stdout.write("\r%.0f %% "%(float(simStep)/float(self.numberOfSimulationSteps)*100.0))
				if int(float(simStep)/float(self.numberOfSimulationSteps)*100.0) == 99:
					print ""

			self.randIndices = np.random.randint(low=0, high=self.numberOfElectronTraps+1, size=(25,2))			

			#TODO: this might slow whole system down
			if rareEarthIndex not in self.randIndices:
				self.randIndices = np.vstack((self.randIndices, np.array(rareEarthIndex)))
				np.random.shuffle(self.randIndices)

			self.randIndices = list(map(tuple, self.randIndices))

			for index in self.randIndices:
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

			# record ground/excited state evolution
			self.electronSystems.recordREstate()
			if not simStep % self.saturationAveraging:
				self.electronSystems.groundStateEvolution.record(simStep, self.electronSystems.rareEarthGroundStateCounter)
				self.electronSystems.excitedStateEvolution.record(simStep, self.electronSystems.rareEarthExcitedStateCounter)
				self.electronSystems.resetEvolutionCounters()
			
			#if not simStep % 10:
			self.recordElectronTrapPopulationDistribution(self.electronSystems.population)

		# after last simulation step
		self.savePopulationEvolution()
		self.saveElectronTrapPopulationDistribution()

		#self.visualize(simStep)

	#--------------------------------------------------------------------------
	def handleRecombination(self):
		# now go through the conduction band's collected electrons and
		# find the ones which can recombine to either an electron trap
		# or to a rare earth.
		while self.cb.numberAvailableElectrons:
			self.electronIndices = self.cb.availableElectronIndices
			#print self.electronIndices
			np.random.shuffle(self.electronIndices)
			electronIndex = tuple(self.electronIndices[0])
			self.possibleRecombinationSlots = np.vstack(np.where(self.electronSystems.population == False))
			self.possibleRecombinationSlots = self.possibleRecombinationSlots.T
			np.random.shuffle(self.possibleRecombinationSlots)
			self.possibleRecombinationSlots = self.possibleRecombinationSlots.T
			self.possibleRecombinationSlots = list(map(tuple, self.possibleRecombinationSlots))
			for esIndex in self.possibleRecombinationSlots:
				esPos = self.electronSystems.getPosition(esIndex)				
				electronPosition = self.cb.getElectronPosition(electronIndex)
				#print esPos, electronPosition
				if np.linalg.norm(esPos - electronPosition) <= self.electronTravelRange:
				#if np.abs(esPos - electronPosition) <= self.electronTravelRange:
					self.cb.donateElectron(electronIndex)
					self.electronSystems.recombine(esIndex)
					break

	#--------------------------------------------------------------------------
	def visualize(self, simStep):
		#if np.sum(self.electronSystems.population) - self.vb.donatedElectronCount == 1:
		#	print "Number electrons is valid."
		#else:
		#	print "Number electrons is INVALID."

		# plot electron distribution
		self.visualizer.visualize(self.pumpBeam, self.stedBeam, self.electronSystems, self.electronSystems._rareEarthIndex, simStep)

		# plot VB electron donation saturation
		#self.vb.evolutionDonatedElectrons.plot()

		#self.electronSystems.groundStateEvolution.plot()
		#self.electronSystems.excitedStateEvolution.plot()
		#self.electronSystems.depletionEvolution.plot()

	#--------------------------------------------------------------------------
	def savePopulationEvolution(self):
		plt.plot(self.electronSystems.groundStateEvolution._x, self.electronSystems.groundStateEvolution._y, label="ground")
		plt.plot(self.electronSystems.excitedStateEvolution._x, self.electronSystems.excitedStateEvolution._y, label="excited")
		plt.legend(loc='best')
		plt.title("Rare earth population evolution, REpos=%g"%(self.electronSystems.getPosition(self.electronSystems._rareEarthIndex)))
		plt.savefig("%sgnd_exc_evolution_REidx_%03d_pump_%.2f_sted_%.2f.png"%(self.savePath, self.positionOfRareEarthCenter, self.pumpAmplitude, self.stedAmplitude))
		plt.close()

	#--------------------------------------------------------------------------
	def recordElectronTrapPopulationDistribution(self, distribution):
		self.electronTrapPopulationDistribution += 1.0*distribution

	#--------------------------------------------------------------------------
	def saveElectronTrapPopulationDistribution(self):
		esXPos = self.electronSystems.x

		host = host_subplot(111, axes_class=AA.Axes)
		plt.subplots_adjust(right=0.75)

		par1 = host.twinx()
		par2 = host.twinx()
		offset = 60
		new_fixed_axis = par2.get_grid_helper().new_fixed_axis
		par2.axis["right"] = new_fixed_axis(loc="right",
		                                    axes=par2,
		                                    offset=(offset, 0))
		par2.axis["right"].toggle(all=True)

		host.set_xlabel("trap position [um]")
		host.set_ylabel("Pump intensity")
		par1.set_ylabel("STED intensity")
		par2.set_ylabel("Population distribution")


		a = self.pumpBeam.profile(esXPos)
		b = self.stedBeam.profile(esXPos)
		p1, = host.plot(esXPos/1.0E-6, a, label="Pump intensity")
		p2, = par1.plot(esXPos/1.0E-6, b, label="STED intensity")
		
		distribution = self.electronTrapPopulationDistribution
		distribution[self.electronSystems._rareEarthIndex] = 0
		p3, = par2.plot(esXPos/1.0E-6, distribution, label="Population distribution")

		host.axis["left"].label.set_color(p1.get_color())
		par1.axis["right"].label.set_color(p2.get_color())
		par2.axis["right"].label.set_color(p3.get_color())

		plt.title("Electron trap population distribution, REpos=%g"%(self.electronSystems.getPosition(self.electronSystems._rareEarthIndex)))
		plt.savefig("%spop_dist_REidx_%03d_pump_%.2f_sted_%.2f.png"%(self.savePath, self.positionOfRareEarthCenter, self.pumpAmplitude, self.stedAmplitude))
		plt.close()

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------

