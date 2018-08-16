from GenBank import *
from PickleBGC import *
from BGC_dictionaries import *

class Controller:

	genBank_obj = GenBank()
	pickle_obj = PickleBGC()
	dictionaries = BGC_dictionaries()

	def __init__(self):
		path_pickle = "" # To be defined
		return self

	def defineMode(int mode):

		if(mode==1):
			genBank_obj.storeGenBankFiles()
			for bgc_obj in genBank_obj.bgcList: 
				pickle_obj.storeBGC(bgc_obj)
				dictionaries.SetDictionaries(bgc_obj)

			dictionaries.exportToText()

		if(mode==2):
			for filename in glob.glob(os.path.join(path, '*.pickle')):
				with open(filename, 'rb') as f:
					bgc = pickle.load(f)
					genBank.bgcList.append(bgc)

		if(mode==3):
			print('NMF')

		if(mode==4):
			print('LDA original')



