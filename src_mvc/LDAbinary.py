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
		self.phi = np.full((g, self.kapa), 0.0001) # Initialise the whole matrix to a small value in order to avoid complications.

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
			row = dictionaries.geneDict[gene]

			if gene in bgc.genes:
				partC += np.sum(bgc.phi[gene]*np.log(self.vita[row]))
			else:
				factor_phi = scipy.special.digamma(np.sum(bgc.gamma))
				array = np.exp(scipy.special.digamma(bgc.gamma) - factor_phi)
				partC += np.sum(array*np.log(1 - self.vita[row]))


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

		####### the algorithm #########
		for bgc in bgcList:
			bgc.gamma_pre = copy.deepcopy(bgc.gamma)

			factor_phi = scipy.special.digamma(np.sum(bgc.gamma))
			array = np.exp(scipy.special.digamma(bgc.gamma) - factor_phi) # this is the exponential gamma factor of the estep
			for gene in dictionaries.geneDict.keys():
				row = dictionaries.geneDict[gene]
				if gene in bgc.genes:
					self.phi[row] = self.vita[row]*array/np.sum(self.vita[row]*array)
					bgc.phi[gene] = copy.deepcopy(self.phi[row])
				else:
					self.phi[row] = (1 - self.vita[row])*array/np.sum((1 - self.vita[row])*array)

				temp += self.phi[row]

			bgc.gamma = self.alpha + temp

			bgc.setPartA(self.alpha)
			bgc.binary_setPartBE(dictionaries.geneDict, self.vita)
			bgc.setPartD()
	##################################
	##### M Step #####################
	##################################

	def MStep(self, bgcList, dictionaries):
		# for each gene from the whole vocabulary
		self.vita_pre = copy.deepcopy(self.vita)

		for gene in dictionaries.geneDict.keys():
			row = dictionaries.geneDict[gene]

			numerator = np.zeros((1,self.kapa), dtype=float)
			denominator = np.zeros((1,self.kapa), dtype=float)

			for bgc in bgcList:
				if gene in bgc.genes:
					numerator += bgc.phi[gene]
					denominator += bgc.phi[gene]
				else:
					factor_phi = scipy.special.digamma(np.sum(bgc.gamma))
					array = np.exp(scipy.special.digamma(bgc.gamma) - factor_phi)
					denominator += (1 - self.vita_pre[row])*array/np.sum((1 - self.vita_pre[row])*array)
			self.vita[row] = (numerator)/(0.000001+denominator)


	def normaliseVita(self):
		# @ Void method
		# for i in range(self.vita.shape[0]):
		# 	# self.vita[i] = self.vita[i]/np.linalg.norm(self.vita[i])
		# 	self.vita[i] = self.vita[i]/np.sum(self.vita[i])
		hello = 'World'