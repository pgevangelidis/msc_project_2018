from GenBank import *
from PickleBGC import *
from BGC_Dictionary import *
from LDA_iterator import *
from LDAmodel import *
from LDAbinary import *
from NMF import *
from MixtureModel import *
from dirCheck import *
from pickle_test import *
from plots_v2 import *
from binaryDatasetCheck import *

class Controller:

	def __init__(self):
		self.exitCond = False
		# self.path_pickle = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\pickle'
		directory = os.path.dirname(os.path.realpath(__file__))
		self.path_pickle = directory + "/pickle/"
		# self.path_pickle = r'/Users/pavlos/Documents/personal/msc_project_2018/src_mvc/pickle/' # To be defined
		folder = dirCheck()
		folder.checkDir(self.path_pickle)
		# self.ldaBGCList = []
		self.genBank_obj = GenBank()
		self.pickle_obj = PickleBGC()
		self.dictionaries = BGC_Dictionary()
		self.health_safety_check = pickle_test()
		self.plotCurve = BGCplot()
		self.nmf = NMF_model()
		self.lda = LDA_model(100) # the number 100 is just an initial number
		self.ldaBinary = LDA_model_binary(100) # the number 100 is just a number
		self.mixModel = MixtureModel(self.dictionaries, 50, 0.01)
		self.modelList = [] # This list should contain 3 models. 1) LDA original 2) LDA binary 3) Mixture Model. Its purpose is to load the plots.

	def defineMode(self, mode):

		if(mode==1):
			print('The files are being processed.\nPlease wait...\n')
			self.genBank_obj.storeGenBankFiles()
			for bgc_obj in self.genBank_obj.bgcList: 
				self.pickle_obj.storeBGC(bgc_obj, self.path_pickle)
				self.dictionaries.setDictionaries(bgc_obj)

			self.dictionaries.exportToText()

		if(mode==2):
			print('The bgc objects are loading.\nPlease wait...\n')
			for filename in glob.glob(os.path.join(self.path_pickle, '*.pickle')):
				bgc_obj = self.pickle_obj.loadBGC(filename)
				# first the bgc will pass the test:
				self.health_safety_check.pickle_loadCheck(bgc_obj)
				self.health_safety_check.testOutput(bgc_obj)
				# Now we are sure it is not empty
				self.genBank_obj.bgcList.append(bgc_obj)
				self.dictionaries.setDictionaries(bgc_obj)
			print('All set. Ready to proceed.\n')
			print("bgc: {}".format(len(self.genBank_obj.bgcList)))

		if(mode==3):
			self.nmf.modelNMF(self.dictionaries)

		if(mode==4):
			print('\nLDA original\n')
			print('{}'.format(len(self.dictionaries.geneDict)))
			lda_obj = LDA_model(len(self.dictionaries.geneDict))
			lda_iter = LDA_iterator(len(self.dictionaries.geneDict), lda_obj, False)
			lda_iter.iterator(self.genBank_obj.bgcList, self.dictionaries)
			self.lda = lda_iter.lda

		if(mode==5):
			print('\nLDA binary\n')
			# create a class which will check if the bgc objects have the right size of phi
			numGenes = len(self.dictionaries.geneDict)
			print('{}'.format(numGenes))

			lda_obj = LDA_model_binary(len(self.dictionaries.geneDict))
			lda_iter = LDA_iterator(len(self.dictionaries.geneDict), lda_obj, True)
			lda_iter.iterator(self.genBank_obj.bgcList, self.dictionaries)
			self.ldaBinary = lda_iter.lda

		if(mode==6):
			print('\nMixture Model\n')
			self.mixModel = MixtureModel(self.dictionaries, 50, 0.0000000001)
			self.mixModel.setSparse()
			self.mixModel.MM_iterator()

		
		if(mode==7):
			self.plotCurve.createPlot(self.lda, self.genBank_obj.bgcList)


		if(mode==8):
			self.exitCond = True

		if(mode==9):
			binarySet = BinaryCheck(self.genBank_obj.bgcList, self.dictionaries.geneDict)
			binarySet.runCheck()



