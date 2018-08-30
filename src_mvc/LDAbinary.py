#This class is the binary version of LDA
# This class is the editted version of E-M step in LDA iterator.
# This class contains the LDA model
import numpy as np
import scipy.special
from BGC import *
from BGC_Dictionary import *
# from lb_check import *
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
		# self.lbTester = lb_check()
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

		# #code added 29/8/2018
		# self.lbTester.partA.append(partA)
		# self.lbTester.partBE.append(partBE)
		# self.lbTester.partD.append(partD)


		for gene in dictionaries.geneDict.keys():
			row = dictionaries.geneDict[gene]

			if gene in bgc.genes:
				partC += np.sum(bgc.phi[gene]*np.log(self.vita[row]))
			else:
				array = np.exp(scipy.special.digamma(bgc.gamma) - scipy.special.digamma(np.sum(bgc.gamma)))
				partC += np.sum(((1 - self.vita[row])*array/np.sum((1 - self.vita[row])*array))*np.log(1 - self.vita[row]))
		# #added code 28/9/2018
		# self.lbTester.partC.append(partC)

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

		# #added code 28/9/2018
		# self.lbTester.check_curve(self.lbTester.partA, "A")
		# self.lbTester.check_curve(self.lbTester.partBE, "BE")
		# self.lbTester.check_curve(self.lbTester.partC, "C")
		# self.lbTester.check_curve(self.lbTester.partD, "D")



	##################################
	##### E Step #####################
	##################################
	def EStep(self, bgcList, dictionaries):
		temp = np.zeros((1,self.kapa), dtype=float)

		####### the algorithm #########
		for bgc in bgcList:
			print("{}\n".format(bgc.name))
			bgc.gamma_pre = copy.deepcopy(bgc.gamma)

			psiGamma = (scipy.special.digamma(bgc.gamma_pre) - scipy.special.digamma(np.sum(bgc.gamma_pre))) # this is the exponential gamma factor of the estep
			for gene in dictionaries.geneDict.keys():
				row = dictionaries.geneDict[gene]
				if gene in bgc.genes:
					self.phi[row] = (self.vita[row]*np.exp(psiGamma))/np.sum(self.vita[row]*np.exp(psiGamma))
					bgc.phi[gene] = copy.deepcopy(self.phi[row])
				else:
					self.phi[row] = (1 - self.vita[row])*np.exp(psiGamma)/np.sum((1 - self.vita[row])*np.exp(psiGamma))

				temp += self.phi[row]

			bgc.gamma = self.alpha + temp
			print("Enter part A.\n")
			bgc.setPartA(self.alpha)
			print("Enter part BE.\n")
			bgc.binary_setPartBE(dictionaries.geneDict, self.vita)
			print("Enter part D.\n")
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
					psiGamma = (scipy.special.digamma(bgc.gamma_pre) - scipy.special.digamma(np.sum(bgc.gamma_pre)))
					denominator += (1 - self.vita[row])*np.exp(psiGamma)/np.sum((1 - self.vita[row])*np.exp(psiGamma))
			self.vita[row] = (0.0000001 + numerator)/(0.00001+denominator)


	def normaliseVita(self):
		# @ Void method
		# for i in range(self.vita.shape[0]):
		# 	# self.vita[i] = self.vita[i]/np.linalg.norm(self.vita[i])
		# 	self.vita[i] = self.vita[i]/np.sum(self.vita[i])
		hello = 'World'