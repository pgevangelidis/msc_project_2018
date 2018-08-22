#This class is the binary version of LDA
# This class is the editted version of E-M step in LDA iterator.
# This class contains the LDA model
import numpy as np
import scipy.special
from BGC import *
from BGC_Dictionary import *
import copy

class LDA_model_binary:

	def __init__(self, g):
		self.kapa = 50 # This is the assumed number of topics. If it is changed here it MUST be changed in BGC file also.
		self.alpha = np.ones((1,self.kapa), dtype=float)
		self.aplha_pre = self.alpha
		self.vita = np.random.rand(g, self.kapa)
		self.vita_pre = self.vita
		self.totalLBound = 0.0
		self.totalLBound_pre = []
		self.phi = {}

	##################################
	#### Lower Bound #################
	##################################
	def lowerBound(self, bgc, dictionaries):

		bgc.lbound_pre.append(bgc.lbound)
		###################
		# Initialisation
		###################
		partA = 0.0
		partBE = 0.0
		partC = 0.0
		partD = 0.0

		partA = bgc.partA
		partBE = bgc.partBE
		partD = bgc.partD


		for gene in dictionaries.geneDict.keys():
			wd = 0
			row = dictionaries.geneDict[gene]
			array = bgc.phi[gene]
			if gene in bgc.genes:
				wd = 1

			partC = np.sum(array[0]*(wd*np.log(self.vita[row] + 0.0001) + (1-wd)*np.log(1-self.vita[row] + 0.0001)))


		temp = partA+partBE+partC+partD

		if np.isinf(temp):
			print('This {} has inf value'.format(bgc.name))
			if(np.isinf(partA)):
				print('Part A has the problem')
			if(np.isinf(partBE)):
				print('Part BE has the problem')
			if(np.isinf(partC)):
				print('Part C has the problem')
			if(np.isinf(partD)):
				print('Part D has the problem')

		bgc.lbound = temp
		self.totalLBound += temp



	##################################
	##### E Step #####################
	##################################
	def EStep(self, bgcList, dictionaries):
		temp = np.zeros((1,50), dtype=float)
		wd = 0
		# bgc.phi_pre = copy.deepcopy(bgc.phi)
		bgc.gamma_pre = copy.deepcopy(bgc.gamma)

		factor_phi = scipy.special.digamma(np.sum(bgc.gamma))
		####### the algorithm #########
		for gene in dictionaries.geneDict.keys():
			wd = 0
			row = dictionaries.geneDict[gene]
			array = np.full((1,50), 0.0001)
			for bgc in bgcList:
				if gene in bgc.genes:
					wd = 1
				# for i in range(self.kapa):
				factor = scipy.special.digamma(bgc.gamma)
				array = ((np.log(self.vita[row]+0.0001)*wd)+(1-wd)*np.log(1-self.vita[row]+0.0001))*np.exp(factor-factor_phi)
				bgc.phi[gene] = array/(np.sum(array)) # this normalising method satisfies the sum to 1


				temp += bgc.phi[gene]

		# tep = np.reshape(temp, (1,temp.shape[0]))
		bgc.gamma = self.alpha + temp

		bgc.setPartA(self.alpha)
		bgc.setPartBE()
		bgc.setPartD()
	##################################
	##### M Step #####################
	##################################

	def MStep(self, bgcList, dictionaries):
		# for each gene from the whole vocabulary
		for gene in dictionaries.geneDict.keys():
			row = dictionaries.geneDict[gene]
			numerator = np.zeros((1,bgc.kapa), dtype=float)
			denominator = np.zeros((1,bgc.kapa), dtype=float)
			# and for each document in the corpus
			for bgc in bgcList:
				# if the gene from the vocabulary exists in the selected bgc then the wd = 1
				if gene in bgc.genes:
					numerator += bgc.phi[gene]
					denominator += bgc.phi[gene]
				denominator += np.log(1-self.vita[row]+0.0001)*np.exp(scipy.special.digamma(bgc.gamma) - scipy.special.digamma(np.sum(bgc.gamma)))

			self.vita[row] = (numerator)/(denominator)


	def normaliseVita(self):
		# @ Void method
		# for i in range(self.vita.shape[0]):
		# 	# self.vita[i] = self.vita[i]/np.linalg.norm(self.vita[i])
		# 	self.vita[i] = self.vita[i]/np.sum(self.vita[i])
		hello = 'World'