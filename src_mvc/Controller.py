from GenBank import *
from PickleBGC import *
from BGC_Dictionary import *
from LDA_iterator import *
from dirCheck import *
from pickle_test import *

class Controller:

	def __init__(self):
		self.exitCond = False
		self.path_pickle = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/pickle/' # To be defined
		folder = dirCheck()
		folder.checkDir(self.path_pickle)

	def defineMode(self, mode):

		genBank_obj = GenBank()
		pickle_obj = PickleBGC()
		dictionaries = BGC_Dictionary()
		health_safety_check = pickle_test()

		if(mode==1):
			genBank_obj.storeGenBankFiles()
			for bgc_obj in genBank_obj.bgcList: 
				pickle_obj.storeBGC(bgc_obj, self.path_pickle)
				dictionaries.setDictionaries(bgc_obj)

			dictionaries.exportToText()

		if(mode==2):
			for filename in glob.glob(os.path.join(self.path_pickle, '*.pickle')):
				bgc = pickle_obj.loadBGC(filename)
				# first the bgc will pass the test:
				health_safety_check.pickle_loadCheck(bgc)
				health_safety_check.testOutput()
				# Now we are sure it is not empty
				genBank_obj.bgcList.append(bgc)

		if(mode==3):
			print('NMF')

		if(mode==4):
			print('\nLDA original\n')
			lda_iter = LDA_iterator(len(dictionaries.geneDict))
			lda_iter.iterator(genBank_obj.bgcList, dictionaries)

		if(mode==6):
			self.exitCond = True



