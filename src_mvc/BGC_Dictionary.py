# This class will store the bgc names to a dictionary with row as a value
# Also,the single genes in a dictionary
# The bgc with the gene list in a dictionary
# The coordinates list of tuples for the sparse matrix
import glob
import os.path
from dirCheck import *

class BGC_Dictionary:

	def __init__(self):
		self.bgcDict = {}
		self.geneDict = {}
		self.bgcGeneDict = {}
		self.coordinates = []
		self.col = 0
		self.row = 0
		self.path = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/product_files/'
		folder = dirCheck()
		folder.checkDir(self.path)

	def setDictionaries(self, bgc):

		if len(self.bgcDict)==0:
			# add the tuple to the dictionary
			self.bgcDict = {bgc.name: self.row}
			self.bgcGeneDict = {bgc.name: bgc.genes}
			self.row += 1
			#print("When I add new BGC the row is: %d" % row)
			# If the list is not empty, check if the BGC already exists.
		else:
			# count how many times does the BGC exists in the list
			if bgc.name not in self.bgcDict:
				# Add the BGC to the dictionary
				self.bgcDict.update({bgc.name: self.row})
				self.bgcGeneDict.update({bgc.name: bgc.genes})
				self.row += 1
				#print("When I add new BGC the row is: %d" % row)
			#else:
			#	print("Houston we 've got a problem!\nThere is a BGC that exists more than once: %s" % bgc.name)


		for geneName in bgc.genes:
			if len(self.geneDict)==0:
				# obviously add the gene to the list
				self.geneDict = {geneName: self.col}
				self.col += 1
			else:
				#first check if the gene exists in the dictionary
				if geneName not in self.geneDict:
					# Now add the new gene to the dictionary
					self.geneDict.update({geneName: self.col})
					self.col += 1
			#	else:
			#		print("Houston we 've got a problem!\nThere is a gene that exists more than once: %s %s" % (geneName,bgc.name))
			# The next line of code stores the coordinates to a list.
			self.setCoordinates(self.bgcDict[bgc.name], self.geneDict[geneName])

	# This method sets the coordinates for the sparse matrix
	def setCoordinates(self, row, col):
		self.coordinates.append((row,col))

	# This method exports the dictionaries and the coordinate list to text files
	def exportToText(self):

		bgc_dict_file = os.path.join(self.path, 'BGC_dictionary.txt')
		gene_dict_file = os.path.join(self.path, 'Gene_dictionary.txt')
		bgc_gene_dict_file = os.path.join(self.path, 'BGC_gene_dictionary.txt')
		coordinates_file = os.path.join(self.path, 'Coordinates_list.txt')

		bgc_out = open(bgc_dict_file,'w')
		bgc_gene_out = open(bgc_gene_dict_file,'w')
		gene_out = open(gene_dict_file,'w')
		coordinate_out = open(coordinates_file,'w')

		bgc_out.write(str(len(self.bgcDict))+'\n')
		bgc_gene_out.write(str(len(self.bgcGeneDict))+'\n')
		for key in self.bgcDict:
		    bgc_out.write(str(key) + '-->' + str(self.bgcDict[key]) + '\n')
		    bgc_gene_out.write(str(key) + ': ' + str(self.bgcGeneDict[key]) + '\n')

		bgc_out.close()
		bgc_gene_out.close()

		gene_out.write(str(len(self.geneDict))+'\n')
		for key in self.geneDict:
		    gene_out.write(str(key) + '-->' + str(self.geneDict[key]) + '\n')

		gene_out.close()

		coordinate_out.write(str(len(self.coordinates))+'\n')
		for item in self.coordinates:
		    coordinate_out.write(str(item)+'\n')

		coordinate_out.close()


		print('The files have been written successfully.\n')
