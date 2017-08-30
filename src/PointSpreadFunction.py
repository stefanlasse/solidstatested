
from multiprocessing import Manager
from threading import Thread

from Simulator import SolidStateStedSimulator

import pickle


class PointSpreadFunction(Thread):
	#--------------------------------------------------------------------------
	def __init__(self, N, REcoord, ETcoord, pumpAmpl, stedAmpl, laserCoord, cs, eTR, savePath):

		super(PointSpreadFunction, self).__init__()

		self.N = N
		self.REcoord = REcoord
		self.ETcoord = ETcoord
		self.pumpAmpl = pumpAmpl
		self.stedAmpl = stedAmpl
		self.laserCoord = laserCoord
		self.crossSections = cs
		self.electronTravelRange = eTR
		self.savePath = savePath

		self.manager = Manager()
		self.resultContainer = self.manager.Queue(maxsize=0)
		self.processList = list()
		self._result = list()

	#--------------------------------------------------------------------------
	def run(self):
		for laserPosition in self.laserCoord:
			print laserPosition
			sim = SolidStateStedSimulator(nSimSteps=self.N, resultContainer=self.resultContainer)
			sim.setupSimulation(REx=self.REcoord[0], REy=self.REcoord[1],
								ETx=self.ETcoord[:,0], ETy=self.ETcoord[:,1],
								pumpAmpl=self.pumpAmpl, stedAmpl=self.stedAmpl,
								laserXpos=laserPosition[0], laserYpos=laserPosition[1],
								cs=self.crossSections, eTR=self.electronTravelRange)
			sim.start()
			self.processList.append(sim)

		for p in self.processList:
			p.join()

		self.saveResult()

	#--------------------------------------------------------------------------
	def saveResult(self):
		while not self.resultContainer.empty():
			self._result.append(self.resultContainer.get())

		with open("%sPSF_pump_%.3f_sted_%.3f.pys"%(self.savePath, self.pumpAmpl, self.stedAmpl), "wb") as f:
			pickle.dump(self._result, f)

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
