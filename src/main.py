# This is the main file
# the main is responsible to read all the .gbk files
# create BGC objects and store all the necessary info
# create two dictionaries for the sparse matrix
# And initialise the LDA


#### packages/libraries/ etc
from Bio import SeqIO
from BGC import *
from BGC_Dictionary import *
from LDAmodel import *
from storeBGC import *
import glob
import os.path
import time
import copy

bgcList = []

if __name__ == '__main__':
	# this is the source folder
	path = r'C:\Users\user\Documents\msc_thesis_2018\gbk\gbk_files' # There is also a test folder under the name \gbk_file_test
	# this is the destination folder
	save_path = r"C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src\geneList"

	####### Timer counter
	start_time = time.time()
	###################

	numberOfFile = 0 # this counter will tell me to which file I have an error.
	# the filename is each BGC.gbk file
	print('********************\n')
	print('Start procedure')
	print('The program is reading all .gbk files and creates BGC objects.')
	print('********************\n')
	for filename in glob.glob(os.path.join(path, '*.gbk')):
		i = 0
		# the geneList is the list of genes from each BGC.gbk file. This list contains all genes of each gene cluster
		# and it will be resetted after each iteration.
		# Actually I used a list with a queue
		end = len(filename)-4
		start = len(filename)-8
		# The output file name must be dynamic and change with each .gbk file
		output_filename = "BGC_GENES_"+filename[start:end]+".txt"
		# I will create the output file in the geneList folder
		completeName = os.path.join(save_path, output_filename)
		# Now open the file to write in it.
		output_handle = open(completeName,'w')

		for gbk_record in SeqIO.parse(filename, "genbank"):


			# for each BGC I will create one object of BGC
			# and then before reading the next BGC store the object
			# to a list object.  Do it at the end.
			bgc = BGC(gbk_record.name)
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

			bgcList.append(bgc) # This line stores the instant object of bgc to a list.

		#print("The file %d was completed successfully." % numberOfFile)
		numberOfFile +=1

	# ############## Test about phi and phi_pre #########################
	# test1 = bgcList[0]
	# test2 = bgcList[23]
	# test3 = bgcList[1023]
	# print('{}:\nthe phi is:\n{}\nphi_pre is:\n{}'.format(test1.name, test1.phi[0], test1.phi_pre[0]))
	# print('{}:\nthe phi is:\n{}\nphi_pre is:\n{}'.format(test2.name, test2.phi[0], test2.phi_pre[0]))
	# print('{}:\nthe phi is:\n{}\nphi_pre is:\n{}'.format(test3.name, test3.phi[0], test3.phi_pre[0]))
	# ###################################################################

	output_handle.close()
	################################
	# Set the dictionaries of BGC, Gene and BGC/Gene
	# Also set the coordinates of
	################################
	dictionaries = BGC_Dictionary()
	for bgc_obj in bgcList:
		dictionaries.setDictionaries(bgc_obj)

	# Use the exportToText method to store the dictionaries
	dictionaries.exportToText()

	readGBK = time.time() - start_time
	tm = time.gmtime(readGBK)
	print('The program has finished. The time needed: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
	print('********************\n')


	print('Now is the LDA phase.')
	######## Apply the LDA Model ##########
	lda = LDA_model(len(dictionaries.geneDict)) # It is important to refer to the same bgc objects I referred above. And modify the same bgc objects

	iteration = 100 # I have picked 100 iteration hoping the LDA will merge before them
	error = 0.01
	errorFlag = False
	loop = 0

	saveBGC = store_BGC()

	while(errorFlag!=True):

		#### loop for all BGC objects the E-step:
		n=0
		print('****** Computation of E - Step ******')
		st = time.time()
		for bgc_obj in bgcList:
			lda.EStep(bgc_obj, dictionaries)
			n+=1

		em = time.time() - st
		tm = time.gmtime(em)
		print('1. Computation of E - Step is done, loop: {} time: {}:{}:{}\n'.format(loop, tm.tm_hour, tm.tm_min, tm.tm_sec))

		#### loop for all BGC objects the M-step:
		n=0
		lda.vita_pre = copy.deepcopy(lda.vita)
		print('****** Computation of M - Step ******')
		st = time.time()
		for bgc_obj in bgcList:
			lda.MStep(bgc_obj, dictionaries)
			n+=1
		mm = time.time() - st
		tm = time.gmtime(mm)
		print('2. Computation of M - Step is done, loop: {} time: {}:{}:{}\n'.format(loop, tm.tm_hour, tm.tm_min, tm.tm_sec))

		# st = time.time()
		# lda.normaliseVita()
		# mm = time.time() - st
		# tm = time.gmtime(mm)
		# print('3. Normalisation of vita matrix is done, loop: {} time: {}:{}:{}\n'.format(loop, tm.tm_hour, tm.tm_min, tm.tm_sec))

		#### Loop for the Lower Bound!
		n=0
		lda.totalLBound_pre.append(lda.totalLBound)
		lda.totalLBound = 0.0
		print("****** Computation of Lower Bound for each BGC ******")
		st = time.time()
		for bgc_obj in bgcList:
			lda.lowerBound(bgc_obj,dictionaries)
			n+=1

		lb = time.time() - st
		tm = time.gmtime(lb)
		print('4. Computation of Lower Bound is done, loop: {} time: {}:{}:{}\n'.format(loop, tm.tm_hour, tm.tm_min, tm.tm_sec))
		h = 0
		if (len(lda.totalLBound_pre)>=1):
			h = len(lda.totalLBound_pre) - 1

		calcError = np.abs(lda.totalLBound - lda.totalLBound_pre[h])
		print('total LB pre: {} post: {}'.format(lda.totalLBound_pre[h],lda.totalLBound))
		#print('total LB printed in main: {}'.format(lda.totalLBound))

		if (calcError <= error) or (loop == iteration):
			errorFlag = True

		print('calcError: {}\t error: {}\nloop: {}\t iteration: {}\n'.format(calcError, error, loop, iteration))
		loop+=1

	loop -= 1
	print('The convergence needed {} loops'.format(loop))
	print('The program stores the BGC objects to text files.\nPlease wait.')


	for bgc_obj in bgcList:
		saveBGC.saveBGCobject(bgc_obj)
	saveBGC.saveLDA(lda,loop)

	elapse_time = time.time() - start_time
	tm = time.gmtime(elapse_time)
	print('The elapsed time is: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
	print('DONE!')
