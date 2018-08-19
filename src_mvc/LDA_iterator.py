# This class was designed to be any interface of iterating between E step and M step, because it can be implemented both in LDA original and binary
from binaryBGC import *
from BGC_Dictionary import *
from LDAmodel import *
from storeBGC import *
import time
import copy

class LDA_iterator:

	def __init__(self, geneDict, lda_obj, binary):
		self.iteration = 60 # I have picked 100 iteration hoping the LDA will merge before them
		self.error = 0.03
		self.errorFlag = False
		self.loop = 0
		self.step = ""
		self.binary = binary
		self.size = geneDict
		self.lda = lda_obj
		self.path = ""
		self.path_b = ""

	def time_status(self, loop, hour, min, sec, s):
		if(s==1):
			self.step = "E - Step"
		if(s==2):
			self.step = "M - Step"
		if(s==3):
			self.step = "Lower Bound"

		print('{}. Computation of {} is done.\nloop: {}, time: {}:{}:{}\n'.format(s, self.step, loop, hour, min, sec))

	# e-step method
	def runEStep(self, bgcList, dictionaries):
		for bgc_obj in bgcList:
			self.lda.EStep(bgc_obj, dictionaries)
	# run m-step method
	def runMStep(self, bgcList, dictionaries):
		for bgc_obj in bgcList:
			self.lda.MStep(bgc_obj, dictionaries)
	# run m-step for binary model
	def runMStepBinary(self, bgcList, dictionaries):
		self.lda.MStep(bgcList, dictionaries)
	# run lower bound equation
	def runLowerBound(self, bgcList, dictionaries):
		for bgc_obj in bgcList:
			self.lda.lowerBound(bgc_obj,dictionaries)
	# choose the right path for each model
	def choosePath(self, binary):
		if binary ==True:
			# windows path
			self.path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\bgc_objects_binary'
			self.path_b = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\lda_object_binary'
			# mac path
			# self.path = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/bgc_objects_binary/')
			# self.path_b = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/lda_object_binary/')
		else:
			# windows path
			self.path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\bgc_objects'
			self.path_b = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\lda_object'
			# mac path
			# self.path = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/bgc_objects/')
			# self.path_b = r'/Users/pavlos/Documents/personal/msc_project/msc_project_2018/src_mvc/lda_object/')


	def iterator(self, bgcList, dictionaries):
		####### Timer counter
		start_time = time.time()
		###################

		while(self.errorFlag!=True):
			#### loop for all BGC objects the E-step:

			print('****** Computation of E - Step ******')
			st = time.time()
			self.runEStep(bgcList, dictionaries)

			em = time.time() - st
			tm = time.gmtime(em)
			self.time_status(self.loop, tm.tm_hour, tm.tm_min, tm.tm_sec, 1)

			#### loop for all BGC objects the M-step:

			self.lda.vita_pre = copy.deepcopy(self.lda.vita)
			
			print('****** Computation of M - Step ******')
			st = time.time()

			if self.binary==True:
				self.runMStepBinary(bgcList, dictionaries)
			else:
				self.runMStep(bgcList, dictionaries)
				self.lda.normaliseVita()

			mm = time.time() - st
			tm = time.gmtime(mm)
			self.time_status(self.loop, tm.tm_hour, tm.tm_min, tm.tm_sec, 2)

			#### Loop for the Lower Bound!

			self.lda.totalLBound_pre.append(self.lda.totalLBound)
			self.lda.totalLBound = 0.0
			print("****** Computation of Lower Bound for each BGC ******")
			st = time.time()
			
			self.runLowerBound(bgcList, dictionaries)

			lb = time.time() - st
			tm = time.gmtime(lb)
			self.time_status(self.loop, tm.tm_hour, tm.tm_min, tm.tm_sec, 3)

			h = 0
			if (len(self.lda.totalLBound_pre)>=1):
				h = len(self.lda.totalLBound_pre) - 1

			calcError = np.abs(self.lda.totalLBound - self.lda.totalLBound_pre[h])
			print('total LB pre: {} post: {}'.format(self.lda.totalLBound_pre[h],self.lda.totalLBound))
			#print('total LB printed in main: {}'.format(lda.totalLBound))

			if (calcError <= self.error) or (self.loop == self.iteration):
				self.errorFlag = True

			print('calcError: {}\t error: {}\nloop: {}\t iteration: {}\n'.format(calcError, self.error, self.loop, self.iteration))
			self.loop+=1

		self.loop -= 1
		print('The convergence needed {} loops'.format(self.loop))
		print('The program stores the BGC objects to text files.\nPlease wait.')

		self.choosePath(self.binary)
		saveBGC = store_BGC(self.path, self.path_b)
		for bgc_obj in bgcList:
			saveBGC.saveBGCobject(bgc_obj)
		saveBGC.saveLDA(self.lda, self.loop)

		elapse_time = time.time() - start_time
		tm = time.gmtime(elapse_time)
		print('The elapsed time is: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
		print('DONE!\n')
