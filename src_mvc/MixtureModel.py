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
		self.pi_pre = np.full((1,self.kapa), (1.0/self.kapa)) # this is the phi parameter 
		self.qiou = np.full((self.ni, self.kapa), (0.0001)) # This is the phi matrix similar to LDA
		self.qiouApprox = np.full((self.ni, self.kapa), 0.0001)
		self.vita_pre = np.zeros((self.di, self.kapa), dtype = float)
		self.vita = np.random.rand(self.di, self.kapa) # This is the equal vita matrix with the gene probabilities
		#This number is more like a safety number if the model cannot converge
		self.iterations = 60
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
				row_gene = self.dictionaries.geneDict[gene] # returns the row of gene
				if gene in self.dictionaries.bgcGeneDict[bgc]:
					array += (np.log(self.vita[row_gene]))
				else:
					array += (np.log(1 - self.vita[row_gene]))
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
		temp = 0.0
		self.pi_pre = copy.deepcopy(self.pi)
		for bgc in self.dictionaries.bgcDict.keys():
			row_bgc = self.dictionaries.bgcDict[bgc]
			log_numerator = np.zeros((1,self.kapa), dtype = float)

			for gene in self.dictionaries.geneDict.keys(): # if the gene exists in the particular bgc then is 1 otherwise 0
				row_gene = self.dictionaries.geneDict[gene] # this line returns the number of row for this gene

				if gene in self.dictionaries.bgcGeneDict[bgc]:
					log_numerator += np.log(self.vita[row_gene] + 0.00001)
				else:
					log_numerator += np.log(1.00001 - self.vita[row_gene])

			# print("numerator: {}\n".format(numerator))
			# print("log_num {}\n".format((log_numerator)))
			# print("num {}\n".format(np.exp(log_numerator)))
			# if np.sum(self.pi_pre*(log_numerator))==0.0:
			# 	self.qiou[row_bgc] = 0.0001
			# 	print('BGC {}\trow number: {}\t the sum is close to zero.\t value: {}'.format(bgc, row_bgc, np.sum((log_numerator))))
			# else:
			self.qiou[row_bgc] = (self.pi_pre*(log_numerator))/np.sum(self.pi_pre*(log_numerator)) # This line normalise the qiou at the same time.
			# self.qiou[row_bgc] = (self.qiou[row_bgc])/np.sum(self.qiou[row_bgc])
			# print("BGC: {}\nqiou:\n{}".format(bgc,self.qiou[row_bgc]))
				
			temp += self.qiou[row_bgc]
		print("show pi before: {}\n".format(self.pi))
		self.pi = temp / self.ni # still a vector. 
		# self.pi = self.pi / np.sum(self.pi) # Normalise to sum to 1.
		print("Show pi: {} N: {}\n".format(self.pi, self.ni))


	def MM_MStep(self):
		# Updating vita
		self.vita_pre = copy.deepcopy(self.vita)

		for gene in self.dictionaries.geneDict.keys():
			row_gene = self.dictionaries.geneDict[gene]
			numerator = np.zeros((1,self.kapa), dtype = float)
			denominator = np.zeros((1,self.kapa), dtype = float)

			for bgc in self.dictionaries.bgcDict.keys():
				row_bgc = self.dictionaries.bgcDict[bgc]

				if gene in self.dictionaries.bgcGeneDict[bgc]:
					numerator += self.qiou[row_bgc]

				denominator += self.qiou[row_bgc]
			self.vita[row_gene] = (numerator)/ (denominator + 0.00001) 

	def MM_iterator(self):
		# This method captures the iterations of EM algorithm which updates the values of the model
		start_time = time.time()
		flag = False
		loop = 0
		# print('vita initial: {}\n'.format(self.vita))
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
			self.calcError = np.square(np.power((self.totalLBound - self.totalLBound_pre[-1]), 2))

			if (self.totalLBound - self.totalLBound_pre[-1]) > 0:
				print("__\n /|\n/")
			else:
				print("  |\n _|_\n \\ /\n  v")

			# sensitive = np.abs(self.calcError - self.calcError_pre)

			if (self.calcError <= self.error) or (loop == self.iterations):
				flag = True
			print('Lower bound:\nPre: {}\tPost:{}\ncalcError: {}\t error: {}\nloop: {}\t iteration: {}\n'.format(self.totalLBound,self.totalLBound_pre[-1], self.calcError, self.error, loop, self.iterations))
			loop += 1
		loop -= 1
		print('The convergence needed {} loops'.format(loop))
		print('The program stores the MM object to text files.\nPlease wait.')

		self.exportMM(loop)
		elapse_time = time.time() - start_time
		tm = time.gmtime(elapse_time)
		print('The elapsed time is: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
		print('DONE!\n')
