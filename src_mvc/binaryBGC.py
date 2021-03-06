#this file will be the old BGC but it will extend the phi matrix to keep record of all genes and probabilites within a dictionary.

# This class is the BGC class
# The class contain information about each BGC
# name, genes, proteins, further useful info and stats

# The class will have a method to store the genes, the definition
# Also to check if there are genes with the same name more than once in the gene cluster
import numpy as np
import scipy.special

class binaryBGC:
	# The constructor will init all the necessary variables

	def __init__(self, n):
		self.name = n
		self.description = ""
		self.source = ""
		self.genes = []
		self.geneDict = {}
		self.doubleGenes = False
		# Those are matrices for the LDA model
		self.kapa = 50 # This is the assumed number of topics. If changed here it MUST change in LDA also.
		self.phi = {} # np.full((1,self.kapa),(1/self.kapa)) # kapa is the number of topics, 1/k
		# self.phi_pre = np.full((1,self.kapa), 0.0)
		self.gamma = np.full((1,self.kapa),((self.kapa + 1)/self.kapa)) # kapa is the number of topics
		self.gamma_pre = np.full((1,self.kapa),0.0) # The gamma = alpha + N/k
		self.number = 0
		self.lbound = 0.0
		self.lbound_pre = []
		self.partA = 0.0
		self.partD = 0.0
		self.partBE = 0.0

	#This method append the gene to the gene Dictionary and the gene list
	def addGene(self, g):
		##### Here I will check if the gene name has space and more than one words.
		##### The purpose is to replace the spaces with dashes

		test = g.split()
		if len(test)>1:
			g = ''
			for tok in test:
				g += str(tok)+'_'

		if g in self.geneDict.keys():
			self.geneDict[g] += 1 # increase the number of gene that exists in a BGC
			self.doubleGenes = True
		else:
			self.geneDict.update({g:1})
			self.genes.append(g) # increase the list with the new gene
			self.number += 1
			# This is the initialisation of phi and gamma part 2
			self.phi.update({g : np.full((1,self.kapa),(1/self.kapa))})
			self.gamma += 1/self.kapa

	# This method checks if there are dublicate genes in
	# the same BGC and returns a list of the gene and the number it exists
	def isDoubleGenes(self):
		dGenes = []
		for key in self.geneDict.keys():
			if self.geneDict[key] > 1:
				dGenes.append((key, self.geneDict[key]))

		return dGenes

	def setPartA(self, alpha):
		self.partA = 0.0
		self.partA = scipy.special.gammaln(np.sum(alpha))
		for i in range(alpha.shape[1]):
			self.partA -= scipy.special.gammaln(alpha[0][i])
			self.partA += (alpha[0][i] - 1)*(scipy.special.digamma(self.gamma[0][i]) - scipy.special.digamma(np.sum(self.gamma)))


	def setPartBE(self):
		self.partBE = 0.0
		for i in self.genes:
			array = self.phi[i]
			for j in range(array.shape[1]):
				self.partBE += array[0][j]*(scipy.special.digamma(self.gamma[0][j]) - scipy.special.digamma(np.sum(self.gamma)))
				self.partBE -= array[0][j]*(np.log(array[0][j]))

	#This method is only for the binary LDA. Because iterates in all genes of the Vocabulary
	# This is also means that the algorithm will be singificantly slower.
	def binary_setPartBE(self, geneDict, vita):
		self.partBE = 0.0 
		psiGamma = scipy.special.digamma(self.gamma) - scipy.special.digamma(np.sum(self.gamma))
		for gene in geneDict:
			row = geneDict[gene]
			if gene in self.genes:
				self.partBE += np.sum(self.phi[gene]*psiGamma)
				self.partBE -= np.sum(self.phi[gene]*(np.log(self.phi[gene])))
			else:
				tempPhi = (1 - vita[row])*np.exp(psiGamma)/np.sum((1 - vita[row])*np.exp(psiGamma))
				self.partBE += np.sum((tempPhi)*(psiGamma))
				self.partBE -= np.sum((tempPhi)*np.log(tempPhi))

	def setPartD(self):
		self.partD = 0.0
		self.partD = (-1)*scipy.special.gammaln(np.sum(self.gamma))
		constant = scipy.special.digamma(np.sum(self.gamma))
		for i in range(self.gamma.shape[1]):
			self.partD += scipy.special.gammaln(self.gamma[0][i])
			self.partD -= (self.gamma[0][i] - 1)*(scipy.special.digamma(self.gamma[0][i]) - constant)

