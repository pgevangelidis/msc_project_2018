import pickle
from BGC import *

class PickleBGC:

	def __init__(self):
		return self

	def storeBGC(self, bgc):
		with open(bgc.name+'.pickle', 'wb') as f:
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

	def loadBGC(self, bgcName):
		try:
			with open( bgcName, 'rb') as f:
				return pickle.load(f)
		except FileNotFoundError:
			print('This name does not exist.')