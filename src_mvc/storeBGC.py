# This file will store to txt file all the information from each BGC
from BGC import *
from LDAmodel import *
import glob
import os.path

class store_BGC:
	def __init__(self, path_a, path_b):
		self.counter = 0
		self.filename = 'Object_'
		self.path = path_a
		self.path_b = path_b
		folder = dirCheck()
		folder.checkDir(self.path)
		folder.checkDir(self.path_b)

	# This is the only method to store all bgc objects to text files
	def saveBGCobject(self, bgc):
		outputName = self.filename+bgc.name+".txt"
		output_file = os.path.join(self.path, outputName)
		output_handler = open(output_file, 'w')
		output_handler.write('name: {}\ndescription: {}\nsource: {}\n\tgenes:\n'.format(bgc.name, bgc.description, bgc.source))
		output_handler.write('\t{}\n'.format(len(bgc.genes)))
		for key in bgc.geneDict.keys():
			output_handler.write('\t{} : {}\n'.format(key, bgc.geneDict[key]))
		#######################################
		#### LDA Variables ####
		output_handler.write('LDA variables\n')
		output_handler.write('kapa: {}\n'.format(bgc.kapa))

		output_handler.write('phi:\n')

		for i in bgc.phi.keys():
			array = bgc.phi[i]
			for j in range(array.shape[0]):
				output_handler.write('\t{}'.format(array[j]))
			output_handler.write('\n')

		output_handler.write('gamma:\n')
		
		for i in range(bgc.gamma.shape[1]):
			output_handler.write('\t{}'.format(bgc.gamma[0][i]))
		output_handler.write('\n')

		output_handler.write('gamma_pre:\n')

		for i in range(bgc.gamma_pre.shape[1]):
			output_handler.write('\t{}'.format(bgc.gamma_pre[0][i]))
		output_handler.write('\n')

		output_handler.write('lower bound:\n{}\n'.format(bgc.lbound))
		output_handler.write('pre lower bound:\n{}\n'.format(bgc.lbound_pre))


		#print('the {} object has been written successfully.'.format(bgc.name))



		output_handler.close()

	def saveLDA(self, lda, loop):
		###################### Total Lower Bound #####################
		outputLB = "totalLowerBound.txt"
		output_file_a = os.path.join(self.path_b, outputLB)
		output_handler_a = open(output_file_a, 'w')
		output_handler_a.write('loop#\tTotal Lower Bound (DESC)\n')
		output_handler_a.write('{}\t{}\n'.format(loop, lda.totalLBound))
		temp = loop
		for i in range(len(lda.totalLBound_pre)-1,0,-1):
			temp -= 1
			output_handler_a.write('{}\t{}\n'.format(temp, lda.totalLBound_pre[i]))

		###################### Vita matrix ##########################
		outputVita = "vita_loop_"+str(loop)+".txt"
		output_file_b = os.path.join(self.path_b, outputVita)
		output_handler_b = open(output_file_b, 'w')
		output_handler_b.write('Gene#\tTopics\n')
		####
		outputVita_pre = "vita_pre_loop_"+str(loop)+".txt"
		output_file_c = os.path.join(self.path_b, outputVita_pre)
		output_handler_c = open(output_file_c, 'w')
		output_handler_c.write('Gene#\tTopics\n')
		####
		for row in range(lda.vita.shape[0]):
			output_handler_b.write("{}\t{}\n".format(row,lda.vita[row]))
			output_handler_c.write("{}\t{}\n".format(row,lda.vita_pre[row]))
