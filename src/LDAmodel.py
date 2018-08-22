# This class contains the LDA model
import numpy as np
import scipy.special
from BGC import *
from BGC_Dictionary import *
import copy

class LDA_model:

	def __init__(self, g):
		self.kapa = 50 # This is the assumed number of topics. If it is changed here it MUST be changed in BGC file also.
		self.alpha = np.ones((1,self.kapa), dtype=float)
		self.aplha_pre = self.alpha
		self.vita = np.random.rand(g, self.kapa)
		self.vita_pre = self.vita
		self.totalLBound = 0.0
		self.totalLBound_pre = []

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

		n = 0
		for gene in bgc.genes:
			row = dictionaries.geneDict[gene]
			for j in range(bgc.phi.shape[1]):
				partC += bgc.phi[n][j]*np.log(self.vita[row][j])
			n += 1

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
	def EStep(self, bgc, dictionaries):
		counter = 0
		bgc.phi_pre = copy.deepcopy(bgc.phi)
		bgc.gamma_pre = copy.deepcopy(bgc.gamma)

		factor_phi = scipy.special.digamma(np.sum(bgc.gamma))
		####### the algorithm #########
		for gene in bgc.genes:
			for i in range(self.kapa):
				row = dictionaries.geneDict[gene]
				factor = scipy.special.digamma(bgc.gamma[0][i])
				bgc.phi[counter][i] = self.vita[row][i]*np.exp(factor-factor_phi)
			bgc.phi[counter] = bgc.phi[counter]/(np.sum(bgc.phi[counter])) # this normalise method satisfies the sum to 1
			counter += 1

		temp = bgc.phi.sum(axis=0)

		tep = np.reshape(temp, (1,temp.shape[0]))
		bgc.gamma = self.alpha + tep

		bgc.setPartA(self.alpha)
		bgc.setPartBE()
		bgc.setPartD()
	##################################
	##### M Step #####################
	##################################

	def MStep(self, bgc, dictionaries):

		# temp = bgc.phi.sum(axis=0)

		n=0 # this counter is for the rows of phi matrix
		for gene in bgc.genes:
			row = dictionaries.geneDict[gene]
			self.vita[row] += bgc.phi[n]
			n += 1

	def normaliseVita(self):

		for i in range(self.vita.shape[0]):
			# self.vita[i] = self.vita[i]/np.linalg.norm(self.vita[i])
			self.vita[i] = self.vita[i]/np.sum(self.vita[i])