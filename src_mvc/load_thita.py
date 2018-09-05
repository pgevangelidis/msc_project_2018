import glob
import os.path
import time
import copy
import re
import numpy as np
from dirCheck import *

if __name__ == "__main__":
	# the directory from where the files will be loaded
	binary_bgc_path = r"C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\bgc_objects_binary"
	bgc_path = r"C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\bgc_objetcs"
	# save path for the cytoscape files
	save_path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\cytoscape_files'
	folder = dirCheck()
	folder.checkDir(save_path)
	############################
	output_style = "styleList_lda.csv"
	s_path = os.path.join(save_path, output_style)
	style_handler = open(s_path, 'w')
	#############################
	output_binary = 'bgc_topic_lda_binary_thita.csv'
	a_path = os.path.join(save_path, output_binary)
	file_handler_a = open(a_path, 'w')

	output_normal = 'bgc_topic_lda_normal_thita.csv'
	b_path = os.path.join(save_path, output_normal)
	file_handler_b = open(b_path, 'w')

	file_handler_a.write('BGC,topic,threshold\n')
	file_handler_b.write('BGC,topic,threshold\n')
	style_handler.write('name,type\n')

	# loop through text files:
	numFiles = 0
	bgcName = ""
	print('Processing...\nPlease wait...\n')
	for filename in glob.glob(os.path.join(binary_bgc_path, '*.txt')):
		try:
			with open(filename, 'r') as f:
				for line in f:
					if 'name:' in line:
						bgcName = line.split()[1]
						style_handler.write('{},BGC\n'.format(bgcName))
					if 'gamma:' in line:
						t = next(f).split()
						t = list(map(float, t))
						gamma = np.asarray(t)
						gamma = gamma/np.sum(gamma)
						for i in range(gamma.shape[0]):
							if (gamma[i]>(np.amax(gamma) - 0.03*np.amax(gamma))):
								file_handler_a.write('{},topic{},{}\n'.format(bgcName,(i+1),gamma[i]))
						break
				numFiles+=1
				
		except:
			print('This name does not exist.')
	# This is necessary for the second directory. The one with the binary bgcs
	numFiles = 0
	bgcName = ""
	for filename in glob.glob(os.path.join(bgc_path, '*.txt')):
		try:
			with open(filename, 'r') as f:
				for line in f:
					if 'name:' in line:
						bgcName = line.split()[1]
						style_handler.write('{},BGC\n'.format(bgcName))
					if 'gamma:' in line:
						t = next(f).split()
						t = list(map(float, t))
						gamma = np.asarray(t)
						gamma = gamma/np.sum(gamma)
						for i in range(gamma.shape[0]):
							if (gamma[i]>(np.amax(gamma) - 0.03*np.amax(gamma))):
								file_handler_b.write('{},topic{},{}\n'.format(bgcName,(i+1),gamma[i]))
						break
				numFiles+=1
				
		except:
			print('This name does not exist.')

	for i in range(50):
		style_handler.write('topic{},topic\n'.format((i+1)))


	print("The Thita process is done.")
	print("Now is turn for the vita matrix.")

	dictionary = {}
	### Open the bgc dictionary and create a style with the genes ###
	directory = os.path.dirname(os.path.realpath(__file__))
	gene_path = directory + "/product_files/Gene_dictionary.txt"
	# try:
	genes = open(gene_path, "r")
	next(genes)
	for line in genes:
		temp_list = re.split(r"(\-->)", line)
		#print(temp_list[0])
		temp_int = re.split(r"(\n)", temp_list.pop())[0]
		#print(int(temp_int))
		dictionary.update({int(temp_int) : temp_list[0]})
		# print(dictionary)
# except:
	# 	print("Dictionary, this name does not exist.")

	binary_bgc_path = directory + "/lda_objects_binary"
	bgc_path = directory + "/lda_objects"

	#############################
	output_binary = 'bgc_topic_lda_binary_vita.csv'
	a_path = os.path.join(save_path, output_binary)
	file_handler_a = open(a_path, 'w')

	output_normal = 'bgc_topic_lda_normal_vita.csv'
	b_path = os.path.join(save_path, output_normal)
	file_handler_b = open(b_path, 'w')

	# loop through text files:
	gene = ""
	print('Processing...\nPlease wait...\n')

	vector_array = np.zeros((1,50))
	vector_list = []
	for filename in glob.glob(os.path.join(binary_bgc_path, 'vita_loop_40.txt')):

		with open(filename, 'r') as f:
			next(f)
			for line in f:
				if "[" in line:
					key = int(line.split()[0])
					gene = dictionary[key]
					# print(gene)
					style_handler.write('{},gene\n'.format(gene))
					# this vector will hold the values of each row, and every time it sees a "[" it will reset.
					vector_array = np.zeros((1,50))
					vector_list = []
					# from this line I am interested for the second part
					temp = re.split(r"[\[\t\r\n\]\s]", line)
					for i in range(len(temp)):
						if temp[i] != "":
							vector_list.append(float(temp[i]))
							# print(temp[i])
				else:
					temp = re.split(r"[\[\t\r\n\]\s]", line)
					for i in range(len(temp)):
						if temp[i] != "":
							vector_list.append(float(temp[i]))
							# print(temp[i])
					if "]" in line:
						temp = re.split(r"[\[\t\r\n\]\s]", line)
						for i in range(len(temp)):
							if temp[i] != "":
								vector_list.append(float(temp[i]))

						vector_array = np.asarray(vector_list)
						vector_array = vector_array/np.sum(vector_array)
						for i in range(vector_array.shape[0]):
							if (vector_array[i]>(np.amax(vector_array))):
								file_handler_a.write('{},topic{},{}\n'.format(gene,(i+1),vector_array[i]))

	# This is necessary for the second directory. The one with the binary bgcs
	for filename in glob.glob(os.path.join(bgc_path, 'vita_loop_60.txt')):
		with open(filename, 'r') as f:
			next(f)
			for line in f:
				if "[" in line:
					key = int(line.split()[0])
					gene = dictionary[key]
					# print(gene)
					style_handler.write('{},gene\n'.format(gene))
					# this vector will hold the values of each row, and every time it sees a "[" it will reset.
					vector_array = np.zeros((1,50))
					vector_list = []
					# from this line I am interested for the second part
					temp = re.split(r"[\[\t\r\n\]\s]", line)
					for i in range(len(temp)):
						if temp[i] != "":
							vector_list.append(float(temp[i]))
							# print(temp[i])
				else:
					temp = re.split(r"[\[\t\r\n\]\s]", line)
					for i in range(len(temp)):
						if temp[i] != "":
							vector_list.append(float(temp[i]))
							# print(temp[i])
					if "]" in line:
						temp = re.split(r"[\[\t\r\n\]\s]", line)
						for i in range(len(temp)):
							if temp[i] != "":
								vector_list.append(float(temp[i]))

						vector_array = np.asarray(vector_list)
						# vector_array = vector_array/np.sum(vector_array)
						for i in range(vector_array.shape[0]):
							vector_array[i]
							if (vector_array[i]>(np.amax(vector_array))):
								file_handler_b.write('{},topic{},{}\n'.format(gene,(i+1),vector_array[i]))
								# print("on.\n")


	style_handler.close()
	file_handler_a.close()
	file_handler_b.close()
	print("All done.")