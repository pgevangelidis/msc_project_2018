# This class reads all genbank files from the folder and writes the necessary attributes
# like genes, name, count genes, duplicates, proteins, 
from LDAmodel import *
from dirCheck import *
import glob
import os.path
from binaryBGC import *
from Bio import SeqIO

class GenBank:


	def __init__(self):
		# # this is the source folder
		# self.path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\gbk_files' # There is also a test folder under the name \gbk_file_test
		# mac path
		directory = os.path.dirname(os.path.realpath(__file__))
		self.path = directory + "\\gbk_test\\"
		# self.path = r'/Users/pavlos/Documents/personal/msc_project_2018/gbk_test/'
		# this is the destination folder
		# self.save_path = r"C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\geneList"
		# # mac path
		directory = os.path.dirname(os.path.realpath(__file__))
		self.save_path = directory + "\\genelist\\"
		# self.save_path = r'/Users/pavlos/Documents/personal/msc_project_2018/src_mvc/genelist/'
		# check the directories
		folder = dirCheck()
		folder.checkDir(self.path)
		folder.checkDir(self.save_path)
		# Declaration of bgcList
		self.bgcList = []


	# If the user select the number 1 option then the controller initiates the 
	# storeGenBankFiles() method to read the GenBank files.
	def storeGenBankFiles(self):
		numberOfFile = 0 # this counter will tell me to which file I have an error.
		# the filename is each BGC.gbk file

		for filename in glob.glob(os.path.join(self.path, '*.gbk')):
			i = 0
			# the geneList is the list of genes from each BGC.gbk file. This list contains all genes of each gene cluster
			# and it will be resetted after each iteration.
			# Actually I used a list with a queue
			end = len(filename)-4
			start = len(filename)-8
			# The output file name must be dynamic and change with each .gbk file
			output_filename = "BGC_GENES_"+filename[start:end]+".txt"
			# I will create the output file in the geneList folder
			completeName = os.path.join(self.save_path, output_filename)
			# Now open the file to write in it.
			output_handle = open(completeName,'w')

			for gbk_record in SeqIO.parse(filename, "genbank"):


				# for each BGC I will create one object of BGC
				# and then before reading the next BGC store the object
				# to a list object.  Do it at the end.
				bgc = binaryBGC(gbk_record.name) # I will change the BGC object to BGCbinary object
				# store the description and the source of the gbk file in to the object
				bgc.description = gbk_record.description
				bgc.source = gbk_record.annotations['source']

				output_handle.write("%s\n" % gbk_record.name)
				flag = False

				for seq_feature in gbk_record.features:
					# There are two types CDS and Gene, both of them contains the gene name. For convinience I used the Gene type.
					if seq_feature.type=='gene':
						flag = True

				for seq_feature in gbk_record.features:
					if flag:
						if seq_feature.type=='gene':
							qualifiers = set(seq_feature.qualifiers.keys())
							if 'gene' in qualifiers:
								bgc.addGene(list(seq_feature.qualifiers['gene'])[0])
								output_handle.write("%s\n" % list(seq_feature.qualifiers['gene'])[0] )
							else:
								key = list(seq_feature.qualifiers.keys())[0]
								bgc.addGene(list(seq_feature.qualifiers[key])[0])
								output_handle.write("%s\n" % list(seq_feature.qualifiers[key])[0] )
					else:
						if seq_feature.type=='CDS':
							qualifiers = set(seq_feature.qualifiers.keys())
							if 'protein_id' in qualifiers:
								bgc.addGene(list(seq_feature.qualifiers['protein_id'])[0])
								output_handle.write("%s\n" % list(seq_feature.qualifiers['protein_id'])[0] )
							else:
								bgc.addGene(list(seq_feature.qualifiers['locus_tag'])[0])
								output_handle.write("%s\n" % list(seq_feature.qualifiers['locus_tag'])[0] )

				self.bgcList.append(bgc) # This line stores the instant object of bgc to a list.

			#print("The file %d was completed successfully." % numberOfFile)
			numberOfFile +=1


			output_handle.close()