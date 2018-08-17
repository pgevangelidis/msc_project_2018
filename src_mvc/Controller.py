from GenBank import *
from PickleBGC import *
from BGC_Dictionary import *
from LDA_iterator import *
from LDAmodel import *
from NMF import *
from dirCheck import *
from pickle_test import *
from plots_v2 import *

class Controller:

	def __init__(self):
		self.exitCond = False
		self.path_pickle = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/pickle/' # To be defined
		folder = dirCheck()
		folder.checkDir(self.path_pickle)
		# self.ldaBGCList = []
		self.genBank_obj = GenBank()
		self.pickle_obj = PickleBGC()
		self.dictionaries = BGC_Dictionary()
		self.health_safety_check = pickle_test()
		self.plotCurve = BGCplot()
		self.nmf = NMF_model()
		self.lda = LDA_model(100)

	def defineMode(self, mode):

		if(mode==1):
			self.genBank_obj.storeGenBankFiles()
			for bgc_obj in self.genBank_obj.bgcList: 
				self.pickle_obj.storeBGC(bgc_obj, self.path_pickle)
				self.dictionaries.setDictionaries(bgc_obj)

			self.dictionaries.exportToText()

		if(mode==2):
			for filename in glob.glob(os.path.join(self.path_pickle, '*.pickle')):
				bgc_obj = self.pickle_obj.loadBGC(filename)
				# first the bgc will pass the test:
				self.health_safety_check.pickle_loadCheck(bgc_obj)
				self.health_safety_check.testOutput()
				# Now we are sure it is not empty
				self.genBank_obj.bgcList.append(bgc_obj)
				self.dictionaries.setDictionaries(bgc_obj)
			# self.ldaBGCList = copy.deepcopy(genBank_obj.bgcList)

		if(mode==3):
			self.nmf.modelNMF(self.dictionaries.coordinates)

		if(mode==4):
			print('\nLDA original\n')
			print('{}'.format(len(self.dictionaries.geneDict)))
			lda_iter = LDA_iterator(len(self.dictionaries.geneDict))
			lda_iter.iterator(self.genBank_obj.bgcList, self.dictionaries)
			self.lda = lda_iter.lda


		
		if(mode==6):
			self.plotCurve.createPlot(self.lda, self.genBank_obj.bgcList)


		if(mode==7):
			self.exitCond = True



