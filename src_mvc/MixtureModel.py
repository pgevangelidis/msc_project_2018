#This class is the Mixture model class
# Uses the sparse matrix from the NMF model
import numpy as np
import copy
import time
from scipy import sparse
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
		self.vita_pre = np.zeros((self.di, self.kapa), dtype = float)
		self.vita = np.random.rand(self.di, self.kapa) # This is the equal vita matrix with the gene probabilities
		#This number is more like a safety number if the model cannot converge
		self.iterations = 3

		directory = os.path.dirname(os.path.realpath(__file__))
		self.path = directory + "/MM_files/"
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
		self.row_ind = []
		self.col_ind = []
		self.mat_values = []
		self.coo = []

    # This method creates the sparse matrix
	def setSparse(self):
		for c in self.dictionaries.coordinates:
			self.row_ind.append(c[0])
			self.col_ind.append(c[1])
			self.mat_values.append(1)
		self.coo = sparse.coo_matrix((self.mat_values, (self.row_ind, self.col_ind))).toarray() # the dimensions are D x V.

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
		for i in range(self.qiou.shape[0]):
			self.partA += np.sum(self.qiou[i]*np.log(self.pi))
		
	def MM_partB(self):
		#calculating part B of Lower Bound
		self.partB = 0.0
		for i in range(self.qiou.shape[0]):
			temp = np.zeros((1,self.kapa))
			for j in range(self.qiou.shape[1]):
				temp[0,j] = np.sum(((self.vita[:,j])**(self.coo[i]))*((1 - self.vita[:,j])**(1 - self.coo[i])))
			temp = self.qiou[i]*temp
			self.partB += np.sum(temp)

	def MM_partC(self):
		#calculating the part C of the lower bound
		self.partC = 0.0
		for i in range(self.qiou.shape[0]):
			self.partC += np.sum(self.qiou[i]*np.log(self.qiou[i]+0.0000001))


	def MM_LBound(self):
		#Calculating the lower bound
		self.totalLBound_pre.append(self.totalLBound)
		self.MM_partA()
		self.MM_partB()
		self.MM_partC() 
		self.totalLBound = self.partA + self.partB - self.partC
 

	def MM_EStep(self):
		# Updating the qiou
		for i in range(self.qiou.shape[0]): # this is actually kapa topics
			for j in range(self.qiou.shape[1]): # this is the BGC corpus
				temp = ((self.vita[:,j])**(self.coo[i]))*((1 - self.vita[:,j])**(1 - self.coo[i]))
				self.qiou[i][j] = self.pi[0,j]*np.prod(temp) # the dimensions are D x k. 
		print("Pre update qiou: {}".format(self.qiou[0]))
		for i in range(self.qiou.shape[0]): # ok, this is the bgc
			self.qiou[i] = (self.qiou[i] + 0.00000001)/(np.sum(self.qiou[i]) + 0.00001) # this line normalise the qiou and also applies the constraint to sum to 1
		print("Post updat: {}".format(self.qiou[0]))


	def MM_MStep(self):
		# Updating vita
		denominator = np.zeros((1,self.kapa))
		for i in range(self.qiou.shape[0]):
			denominator += self.qiou[i]

		numerator = np.zeros((1,self.kapa))
		for i in range(self.vita.shape[0]):
			for j in range(self.qiou.shape[1]):
				numerator[0,j] = np.sum(self.qiou[:,j]*self.coo[:,i])
			self.vita[i] = (numerator + 0.0000001)/(denominator + 0.00001)

		# update pi as well
		print("Pi pre updtaing: {}\n".format(self.pi))
		for i in range(self.qiou.shape[0]):
			self.pi += self.qiou[i]
		self.pi = self.pi/self.di
		print("Pi post: {}\n".format(self.pi))


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
