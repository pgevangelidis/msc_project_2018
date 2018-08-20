#This class is the Mixture model class
# Uses the sparse matrix from the NMF model
import numpy as np
import copy

class MixtureModel:

	def __init__(self, dictionaries, kapa, err):
		#This number wants the model to calibrate. Kind of. 
		self.error = err
		self.ni = len(dictionaries.bgcDict.keys()) 
		self.kapa = kapa
		self.di = len(dictionaries.geneDict.keys())
		self.pi = np.full((1,self.kapa), (1/self.kapa)) # this is the phi parameter 
		self.qiou = np.zeros((self.ni, self.kapa), dtype=float) # This is the phi matrix similar to LDA
		self.vita = np.random.rand(self.di, self.kapa) # This is the equal vita matrix with the gene probabilities
		#This number is more like a safety number if the model cannot converge
		self.iterations = 100
		# Here I store the dictionaries
		self.dictionaries = dictionaries # this object contains 3 important dictionaries 1) bgcDict 2) geneDict 3) bgcGeneDict
		# mac path
        # self.path = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/NMF_files/'
        # windows path
        self.path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\MM_files'
        folder = dirCheck()
        folder.checkDir(self.path)
        # part of the equation
        self.partA = 0.0
        self.partB = 0.0
        self.partC = 0.0
        self.LB = 0.0
        self.LB_pre = 0.0

	def exportMM(self):

        Qiou_file = os.path.join(self.path, 'Qiou_file.txt')
        Beta_file = os.path.join(self.path, 'Beta_file.txt')

        Qiou_out = open(Qiou_file,'w')
        Beta_out = open(Beta_file,'w')

        for i in range(self.qiou.shape[0]):
            Qiou_out.write('{}'.format(self.qiou[i]))

        for i in range(self.vita.shape[0]):
            Beta_out.write('{}'.format(self.vita[i]))

        Qiou_out.close()
        Beta_out.close()

    def MM_partA(self):
    	# calculate part A of Lower Bound
    	self.partA = 0.0
    	for bgc in self.dictionaries.bgcDict.keys():
    		row_bgc = self.dictionaries.bgcDict[bgc]
    		self.partA += np.sum(self.qiou[bgc]*np.log(self.pi))
		
	def MM_partB(self):
		#calculating part B of Lower Bound
		self.partB = 0.0
		for bgc in self.dictionaries.bgcDict.keys():			
			row_bgc = self.dictionaries.bgcDict[bgc] # return the row of bgc
			array = np.zeros((1,self.kapa), dtype = float)
			for gene in self.dictionaries.geneDict.keys():
				wd = 0
				row_gene = self.dictionaries.geneDict[gene] # returns the row of gene
				if gene in self.dictionaries.bgcGeneDict[bgc]:
					wd = 1
				array += (wd*np.log(self.vita[row_gene]) + (1-wd)*np.log(1 - self.vita[row_gene]))
			self.partB += np.sum(self.qiou[row_bgc]*array) # this should return a scalar value. I hope.

	def MM_partC(self):
		#calculating the part C of the lower bound
		self.partC = 0.0
		for bgc in self.dictionaries.bgcDict.keys():
    		row_bgc = self.dictionaries.bgcDict[bgc]
    		self.partC += np.sum(self.qiou[row_bgc]*np.log(self.qiou[row_bgc]))


    def MM_LBound(self):
    	#Calculating the lower bound
    	self.LB_pre = copy.deepcopy(self.LB)
    	self.MM_partA()
    	self.MM_partB()
    	self.MM_partC() 
    	self.LB = self.partA + self.partB - self.partC



	def MM_EStep(self):
		# Updating the qiou
		for bgc in self.dictionaries.bgcDict.keys():
			row_bgc = self.dictionaries.bgcDict[bgc]
			numerator = np.zeros((1,self.kapa), dtype = float)
			for gene in self.dictionaries.geneDict.keys(): # if the gene exists in the particular bgc then is 1 otherwise 0
				wd = 0
				row_gene = self.dictionaries.geneDict[gene] # this line reutrns the number of row for this gene
				if gene in self.dictionaries.bgcGeneDict[bgc]:
					wd = 1
				numerator += wd*np.log(self.vita[row_gene]) + (1-wd)*np.log(1 - self.vita[row_gene])
			self.qiou[row_bgc] = (self.pi * numerator) / np.sum(self.pi * numerator) # This line normalise the qiou at the same time.


		# Updating the pi
		for bgc in self.dictionaries.bgcDict.keys():
			row_bgc = self.dictionaries.bgcDict[bgc]
			self.pi += self.qiou
		self.pi = self.pi / self.ni # still a vector. Normalised too.


	def MM_MStep(self):
		# Updating vita
		for gene in self.dictionaries.geneDict.keys():
			row_gene = self.dictionaries.geneDict[gene]
			numerator = np.zeros((1,self.kapa), dtype = float)
			denominator = np.zeros((1,self.kapa), dtype = float)
			for bgc in self.dictionaries.bgcDict.keys():
				wd = 0
				row_bgc = self.dictionaries.bgcDict[bgc]
				if gene in self.dictionaries.bgcGeneDict[bgc]:
					wd = 1 
				numerator += self.qiou[bgc]*wd
				denominator += self.qiou[bgc]
			self.vita[row_gene] = numerator / denominator  

	def MM_iterator(self):
		# This method captures the iterations of EM algorithm which updates the values of the model
		start_time = time.time()
		flag = False
		loop = 0
		calcError = 0.0
		while(flag!=True):

			self.MM_EStep()
			self.MM_MStep()
			self.MM_LBound()

			calcError = np.abs(self.LB - self.LB_pre)

			if (calcError <= self.error) or (loop == self.iterations):
				flag = True
			print('calcError: {}\t error: {}\nloop: {}\t iteration: {}\n'.format(calcError, self.error, loop, self.iterations))
			loop += 1
		loop -= 1
		print('The convergence needed {} loops'.format(loop))
		print('The program stores the MM object to text files.\nPlease wait.')

		self.exportMM()
		elapse_time = time.time() - start_time
		tm = time.gmtime(elapse_time)
		print('The elapsed time is: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
		print('DONE!\n')
