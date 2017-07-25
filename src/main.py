

from Simulator import SolidStateStedSimulator

s = SolidStateStedSimulator(nSimSteps=5E5, nET=300, eTR=25E-9, posRE=0.0, centerET=0.0, spanET=2E-6)
s.start()

s.join()

s.visualize()