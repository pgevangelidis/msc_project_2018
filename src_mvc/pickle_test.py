# This class exists to check if the pickle files have been loaded correctly

class pickle_test:

	def __init__(self):
		self.topic = 50
		self.flag = True
		self.tests = 0

	def testOutput(self):
		if self.flag == False:
			print('Something went wrong.\nThe load was NOT successfull.\n')
		else:
			print('All set.\n')


	def is_empty(self, any_structure):
		if any_structure:
			# print('Structure is not empty.')
			return False
		else:
			# print('Structure is empty.')
			return True

	def pickle_loadCheck(self, bgc):
		self.tests = 0
		# 6 tests have been counted. If the tests count is not zero then the load of BGC was not successfull
		if not bgc.name:
			self.tests += 1

		if self.is_empty(bgc.genes)==True:
			self.tests += 1

		if self.is_empty(bgc.geneDict)==True:
			self.tests += 1

		if bgc.kapa != self.topic:
			self.tests += 1

		if bgc.gamma.shape[1] != self.topic and self.is_empty(bgc.gamma)==True:
			self.tests += 1

		if bgc.phi.shape[1] != self.topic and self.is_empty(bgc.phi)==True:
			self.tests += 1


		if self.tests > 0:
			self.flag = False


	






