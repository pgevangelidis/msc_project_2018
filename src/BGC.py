# This class is the BGC class
# The class contain information about each BGC
# name, genes, proteins, further useful info and stats

# The class will have a method to store the genes, the definition
# Also to check if there are genes with the same name more than once in the gene cluster
import numpy as np

class BGC:
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
		self.phi = np.full((1,self.kapa),(1/self.kapa)) # kapa is the number of topics, 1/k
		self.phi_pre = np.full((1,self.kapa), 0.0)
		self.gamma = np.full((1,self.kapa),(51/self.kapa)) # kapa is the number of topics
		self.gamma_pre = np.full((1,self.kapa),0.0) # The gamma = alpha + N/k
		self.number = 0
		self.lbound = 0.0
		self.lbound_pre = []

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
			if (len(self.geneDict) > 1):
				self.phi = np.vstack((self.phi, np.full((1,self.kapa),(1/self.kapa))))
				self.phi_pre = np.vstack((self.phi_pre, np.full((1,self.kapa), 0.0)))
				self.gamma += 1/self.kapa

	# This method checks if there are dublicate genes in
	# the same BGC and returns a list of the gene and the number it exists
	def isDoubleGenes(self):
		dGenes = []
		for key in self.geneDict.keys():
			if self.geneDict[key] > 1:
				dGenes.append((key, self.geneDict[key]))

		return dGenes
