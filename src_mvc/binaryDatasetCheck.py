### This script will check if the data set is basically binary or not.
import glob
import os.path
import time
import copy
import re
import numpy as np
import pickle
from binaryBGC import *
from PickleBGC import *
from BGC_Dictionary import *

class BinaryCheck:

	def __init__(self, list, gDict):
		self.doubleList = []
		self.multipleList = []
		self.bgcList = list
		self.geneDict = gDict

	def runCheckForDoubles(self):
		for bgc in self.bgcList:
			dGenes = bgc.isDoubleGenes()
			if dGenes:
				for item in dGenes:
					self.doubleList.append((bgc.name, item[0], item[1]))

	def runCheckForMultiples(self):
		print(len(self.bgcList))
		for gene in self.geneDict.keys():
			counter = 0
			tempList = []
			for bgc in self.bgcList:
				for g in bgc.genes:
					if gene==g:
						# print("gene: {} g: {}".format(gene, g))
						counter += 1
						tempList.append(bgc.name)
						break
			if counter>1:
				# print("gene: {}, list: {}, counter: {}".format(gene, tempList, counter))		
				self.multipleList.append((gene, tempList, counter))

	def runAnalysis(self):
		total_genes = len(self.geneDict)
		multiple_genes = len(self.multipleList)
		double_genes = len(self.doubleList)
		print("The total number of genes are: {}".format(total_genes)) 
		print("The genes that exist in more than one BGC are: {} {}%".format(multiple_genes, round(multiple_genes*100.0/total_genes, 2)))
		print("The genes that exist more than once in a single BGC are: {} {}%".format(double_genes, round(double_genes*100.0/total_genes, 2)))

		for item in self.multipleList:
			counter = 0
			for sec_item in self.doubleList:
				# print("{} == {}".format(item[0], sec_item[1]))
				if item[0]==sec_item[1]:
					counter += 1

		print("The union of multiple and duplicate genes: {}".format(round(counter*100.0/double_genes, 2)))


	def runCheck(self):
		self.runCheckForDoubles()
		self.runCheckForMultiples()
		self.runAnalysis()


