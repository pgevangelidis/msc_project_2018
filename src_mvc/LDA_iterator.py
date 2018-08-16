# This class was designed to be any interface of iterating between E step and M step, because it can be implemented both in LDA original and binary
from BGC import *
from BGC_Dictionary import *
from LDAmodel import *
from storeBGC import *
import time
import copy

class LDA_iterator:

	def __init__(self, geneDict):
		self.iteration = 10 # I have picked 100 iteration hoping the LDA will merge before them
		self.error = 0.01
		self.errorFlag = False
		self.loop = 0
		self.step = ""
		self.size = geneDict

	def time_status(self, loop, hour, min, sec, s):
		if(s==1):
			self.step = "E - Step"
		if(s==2):
			self.step = "M - Step"
		if(s==3):
			self.step = "Lower Bound"

		print('{}. Computation of {} is done.\nloop: {}, time: {}:{}:{}\n'.format(s, self.step, loop, hour, min, sec))


	def iterator(self, bgcList, dictionaries):
		lda = LDA_model(self.size)
		saveBGC = store_BGC()

		####### Timer counter
		start_time = time.time()
		###################

		while(self.errorFlag!=True):
			#### loop for all BGC objects the E-step:
			n=0
			print('****** Computation of E - Step ******')
			st = time.time()
			for bgc_obj in bgcList:
				print(bgc_obj.names)
				lda.EStep(bgc_obj, dictionaries)
				n+=1

			em = time.time() - st
			tm = time.gmtime(em)
			self.time_status(self.loop, tm.tm_hour, tm.tm_min, tm.tm_sec, 1)


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
			self.time_status(self.loop, tm.tm_hour, tm.tm_min, tm.tm_sec, 2)

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
			self.time_status(self.loop, tm.tm_hour, tm.tm_min, tm.tm_sec, 3)

			h = 0
			if (len(lda.totalLBound_pre)>=1):
				h = len(lda.totalLBound_pre) - 1

			calcError = np.abs(lda.totalLBound - lda.totalLBound_pre[h])
			print('total LB pre: {} post: {}'.format(lda.totalLBound_pre[h],lda.totalLBound))
			#print('total LB printed in main: {}'.format(lda.totalLBound))

			if (calcError <= self.error) or (self.loop == self.iteration):
				self.errorFlag = True

			print('calcError: {}\t error: {}\nloop: {}\t iteration: {}\n'.format(calcError, self.error, self.loop, self.iteration))
			self.loop+=1

		self.loop -= 1
		print('The convergence needed {} loops'.format(self.loop))
		print('The program stores the BGC objects to text files.\nPlease wait.')

		for bgc_obj in bgcList:
			saveBGC.saveBGCobject(bgc_obj)
		saveBGC.saveLDA(lda, self.loop)

		elapse_time = time.time() - start_time
		tm = time.gmtime(elapse_time)
		print('The elapsed time is: {}:{}:{}\n'.format(tm.tm_hour, tm.tm_min, tm.tm_sec))
		print('DONE!\n')