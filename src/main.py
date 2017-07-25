
from ElectronicSystems import ElectronTrap, RareEarth
from Crystal import ConductionBand, ValenceBand
from LaserProfiles import PumpBeam, StedBeam
from Visualizer import Visualizer

import numpy as np
import matplotlib.pyplot as plt
import random
random.seed()
import copy


# configuration part
numberOfSimulationSteps = int(2E6)
numberOfElectronTraps = 300
electronTravelRange = 25E-9
electronTrapXPosition = np.linspace(-1E-6, 1E-6, numberOfElectronTraps)
positionOfRareEarthCenter = 0.0

# set up pump and sted beam
pumpBeam = PumpBeam(x=0.0, amplitude=0.005, wavelength=470E-9, numAperture=1.3)
stedBeam = StedBeam(x=0.0, amplitude=0.01,  wavelength=600E-9, numAperture=1.3)

# set up conduction and valence band
cb = ConductionBand()
vb = ValenceBand()

# build up electron traps
electronTraps = [ElectronTrap(x=electronTrapXPosition[i]) for i in range(numberOfElectronTraps)]

# set up rare earth
rareEarth = RareEarth(x=positionOfRareEarthCenter)

# calculate laser intensities at each electron trap and rare earth position
for et in electronTraps:
	et.localPumpPower = pumpBeam.profile(et.x)
	et.localStedPower = stedBeam.profile(et.x)

rareEarth.localPumpPower = pumpBeam.profile(rareEarth.x)
rareEarth.localStedPower = stedBeam.profile(rareEarth.x)

# set up transition probabilities for the electron traps
for et in electronTraps:
	et.probabilityPumpIonize = 0.5*pumpBeam.profile(et.x)
	et.probabilityStedIonize = et.probabilityPumpIonize + 0.5*stedBeam.profile(et.x)
	#et.probabilityRecombine  = 0.0

# set up transition probabilities for the rare earth
rareEarth.probabilityDecay       = 0.01
rareEarth.probabilityPumpIonize  = rareEarth.probabilityDecay       + pumpBeam.profile(rareEarth.x)
rareEarth.probabilityStedIonize  = rareEarth.probabilityPumpIonize  + stedBeam.profile(rareEarth.x)
rareEarth.probabilityPumpExcite  = rareEarth.probabilityStedIonize  + 0.2*pumpBeam.profile(rareEarth.x)
rareEarth.probabilityPumpRestore = rareEarth.probabilityPumpExcite  + pumpBeam.profile(rareEarth.x)
rareEarth.probabilityStedDeplete = rareEarth.probabilityPumpRestore + 0.2*stedBeam.profile(rareEarth.x)
#rareEarth.probabilityRecombine   = 0.0

electronTraps.append(copy.deepcopy(rareEarth))
del rareEarth


# start simulation loop
for simStep in range(numberOfSimulationSteps):
	if simStep % 100000 == 0:
		random.seed()
		print "\r%.1f %%"%(float(simStep)/float(numberOfSimulationSteps)*100.0)

	random.shuffle(electronTraps)

	for element in electronTraps:
		randomNumber = random.random()

		if type(element) is RareEarth:
			if randomNumber <= element.probabilityDecay:
				element.decay()

			elif randomNumber > element.probabilityDecay and randomNumber <= element.probabilityPumpIonize:
				cb.absorbElectron(element.ionize())

			elif randomNumber > element.probabilityPumpIonize and randomNumber <= element.probabilityStedIonize:
				cb.absorbElectron(element.ionize())

			elif randomNumber > element.probabilityStedIonize and randomNumber <= element.probabilityPumpExcite:
				element.excite()

			elif randomNumber > element.probabilityPumpExcite and randomNumber <= element.probabilityPumpRestore:
				if not element.isPopulated:
					element.restore(vb.donateElectron(element.x, electronTravelRange))
					vb.recordEvolution(simStep)

			elif randomNumber > element.probabilityPumpRestore and randomNumber <= element.probabilityStedDeplete:
				element.deplete()

			else:
				pass
				

		elif type(element) is ElectronTrap:
			if randomNumber <= element.probabilityStedIonize:
				cb.absorbElectron(element.ionize())

			else:
				pass

		else:
			pass

	#==========================================================================
	# recombination part
	#==========================================================================

	# now go through the conduction band's collected electrons and
	# find the ones which can recombine to either an electron trap
	# or to a rare earth.
	cb.shuffleElectrons()
	while cb.availableElectrons != 0:
		# get electron from shuffled list
		e = cb.getElectron(0)

		# search for an electron trap/RE which is within the electrons
		# travel range and where it can recombine if the trap/RE is currently ionized
		for element in electronTraps:
			if np.abs(element.x - e.x) <= e.travelRange and not element.isPopulated:
				element.recombine(cb.donateElectron(0))
				break

# sort electronic systems by x-position
electronTraps.sort(key=lambda element: element.x, reverse=True)


#for element in electronTraps:
#	print str(element), "is populated:", element.isPopulated

v = Visualizer()
v.visualize(pumpBeam, stedBeam, electronTraps)

vb.plotElectronEvolution()

