#This class is the Mixture model class
# Uses the sparse matrix from the NMF model
import numpy as np
import copy
import time
import os
from dirCheck import *
from storeBGC import *

class MixtureModel:

	def __init__(self, diction, k, err):
		# Here I store the dictionaries
		self.dictionaries = diction # this object contains 3 important dictionaries 1) bgcDict 2) geneDict 3) bgcGeneDict
		#This number wants the model to calibrate. Kind of. 
		self.error = err
		self.ni = len(diction.bgcDict.keys()) 
		self.kapa = k
		self.di = len(diction.geneDict.keys())
		self.pi = np.full((1,self.kapa), (1.0/self.kapa)) # this is the phi parameter 
		self.qiou = np.full((self.ni, self.kapa), (0.0001)) # This is the phi matrix similar to LDA
		self.qiouApprox = np.full((self.ni, self.kapa), 0.0001)
		self.vita_pre = np.zeros((self.di, self.kapa), dtype = float)
		self.vita = np.random.rand(self.di, self.kapa) # This is the equal vita matrix with the gene probabilities
		#This number is more like a safety number if the model cannot converge
		self.iterations = 25
		# mac path
		directory = os.path.dirname(os.path.realpath(__file__))
		self.path = directory + "/MM_files/"
		# self.path = r'/Users/pavlos/Documents/personal/msc_project_2018/src_mvc/MM_files/'
		# windows path
		# self.path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\MM_files'
		folder = dirCheck()
		folder.checkDir(self.path)
		# part of the equation
		self.partA = 0.0
		self.partB = 0.0
		self.partC = 0.0
		self.totalLBound = 0.0
		self.totalLBound_pre = []
		self.calcError = 0.0
		self.calcError_pre = 0.0
		self.sensitive = 0.01

	def exportMM(self, loop):

		saveMixture = store_BGC(self.path, self.path)
		saveMixture.saveLDA(self, loop)

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
			self.partA += np.sum(self.qiou[row_bgc]*np.log(self.pi))
		
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
		self.totalLBound_pre.append(self.totalLBound)
		self.MM_partA()
		self.MM_partB()
		self.MM_partC() 
		self.totalLBound = self.partA + self.partB - self.partC
 

	def MM_EStep(self):
		# Updating the qiou
		for bgc in self.dictionaries.bgcDict.keys():
			row_bgc = self.dictionaries.bgcDict[bgc]
			numerator = np.zeros((1,self.kapa), dtype = float)
			# print("pre\n{}\n".format(numerator))
			first = 0
			for gene in self.dictionaries.geneDict.keys(): # if the gene exists in the particular bgc then is 1 otherwise 0
				row_gene = self.dictionaries.geneDict[gene] # this line returns the number of row for this gene
				# print(self.vita[row_gene])
				if gene in self.dictionaries.bgcGeneDict[bgc]:
					if first==0:
						numerator = 1*(self.vita[row_gene])
					else:
						numerator = numerator*(self.vita[row_gene])
				else:
					if first==0:
						numerator = 1*(1.0001 - self.vita[row_gene])
					else:
						numerator = numerator*(1.0001 - self.vita[row_gene])
				# print("post\n{}\n".format(numerator)) # The smoothness is for the divide by zero
			# print("numerator: {}".format(numerator))
			# print("pi {}".format(self.pi))
			self.qiou[row_bgc] = (self.pi * numerator) / np.sum(self.pi * numerator) # This line normalise the qiou at the same time.
			# print("BGC: {}\nqiou:\n{}".format(bgc,self.qiou[row_bgc]))
				
			

		# Updating the pi
		for bgc in self.dictionaries.bgcDict.keys():
			row_bgc = self.dictionaries.bgcDict[bgc]
			self.pi += self.qiou[row_bgc]
		self.pi = self.pi / self.ni # still a vector. Normalised too.


	def MM_MStep(self):
		# Updating vita
		self.vita_pre = copy.deepcopy(self.vita)

		for gene in self.dictionaries.geneDict.keys():
			row_gene = self.dictionaries.geneDict[gene]
			numerator = np.zeros((1,self.kapa), dtype = float)
			denominator = np.zeros((1,self.kapa), dtype = float)
			for bgc in self.dictionaries.bgcDict.keys():
				wd = 0
				row_bgc = self.dictionaries.bgcDict[bgc]
				if gene in self.dictionaries.bgcGeneDict[bgc]:
					wd = 1 
				numerator += self.qiou[row_bgc]*wd
				denominator += self.qiou[row_bgc]
			self.vita[row_gene] = (numerator + 0.0001)/ (denominator + 0.005) 

	def MM_iterator(self):
		# This method captures the iterations of EM algorithm which updates the values of the model
		start_time = time.time()
		flag = False
		loop = 0
		while(flag!=True):
			print('loop {}\n'.format(loop))
			print('\n***** Computation of E - Step *****\n')
			self.MM_EStep()
			print('\n***** Computation of M - Step *****\n')
			print('This will take some time')
			self.MM_MStep()
			print('\n***** Computation of Lower Bound *****\n')
			self.MM_LBound()

			self.calcError_pre = copy.deepcopy(self.calcError)
			self.calcError = np.abs(self.totalLBound - self.totalLBound_pre[-1])
			sensitive = np.abs(self.calcError - self.calcError_pre)

			if ((self.calcError <= self.error) and (sensitive <= self.error)) or (loop == self.iterations):
				flag = True
			print('calcError: {}\t sensitive: {}\t error: {}\nloop: {}\t iteration: {}\n'.format(self.calcError, sensitive, self.error, loop, self.iterations))
			loop += 1
		loop -= 1
		print('The convergence needed {} loops'.format(loop))
		print('The program stores the MM object to text files.\nPlease wait.')

		self.exportMM(loop)
		elapse_time = time.time() - start_time
		tm = time.gmtime(elapse_time)
		print('The elapsed time is: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
		print('DONE!\n')
