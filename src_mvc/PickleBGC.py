import pickle
from BGC import *

class PickleBGC:

	def storeBGC(self, bgc, dir):
		with open(dir + '\\' + bgc.name + '.pickle', 'wb') as f:
			pickle.dump(bgc, f, pickle.HIGHEST_PROTOCOL)

	def loadBGC(self, bgcName):
		try:
			with open(bgcName, 'rb') as f:
				return pickle.load(f)
		except:
			print('This name does not exist.')

	def storeModel(self, model, name, dir):
		with open(dir + "\\" + name + ".pickle", "wb") as f:
			pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)
