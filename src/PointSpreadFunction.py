
from multiprocessing import Queue, TimeoutError
import Queue as queue
from Simulator import SolidStateStedSimulator
from threading import Thread

import time

class PointSpreadFunction(Thread):
	#--------------------------------------------------------------------------
	def __init__(self, N, REcoord, ETcoord, pumpAmpl, stedAmpl, cs, eTR):
		self.N = N
		self.REcoord = REcoord
		self.ETcoord = ETcoord
		self.pumpAmpl = pumpAmpl
		self.stedAmpl = stedAmpl
		self.crossSections = cs
		self.electronTravelRange = eTR

		self.resultContainer = Queue()
		self.queueGrabber    = QueueGrabber(len(REcoord), self.resultContainer )
		self.processList = list()

		Thread.__init__(self)

	#--------------------------------------------------------------------------
	def setup(self):
		for rePos in self.REcoord:
			sim = SolidStateStedSimulator(nSimSteps=self.N, resultContainer=self.resultContainer)
			sim.setupSimulation(REx=rePos[0], REy=rePos[1],
								ETx=self.ETcoord[:,0], ETy=self.ETcoord[:,1],
								pumpAmpl=self.pumpAmpl, stedAmpl=self.stedAmpl,
								cs=self.crossSections, eTR=self.electronTravelRange)

			self.processList.append(sim)

	#--------------------------------------------------------------------------
	def run(self):
		for p in self.processList:
			p.start()

		self.waitUntilFinished()
		return

	#--------------------------------------------------------------------------
	def waitUntilFinished(self):
		self.queueGrabber.start()
		self.queueGrabber.join()
		for p in self.processList:
			print "Trying to join process.", p.pid
			p.join()
			print "Joint process.", p.pid

	#--------------------------------------------------------------------------
	def saveResult(self):
		pass

	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------

class QueueGrabber(Thread):
	def __init__(self, n, q):
		self.n =n
		self.q = q
		self.data = list()
		self.stop = False
		Thread.__init__(self)

	def get(self):
		try:
			data = self.q.get(timeout=0.1)
		except TimeoutError as err:
			pass
		except queue.Empty as err:
			pass			
		else:
			self.data.append(data)			

	def stop(self):
		self.stop = True

	def run(self):
		while (self.q.qsize() < self.n) and not self.stop:
			self.get()
			time.sleep(0.5)

		self.stop = False