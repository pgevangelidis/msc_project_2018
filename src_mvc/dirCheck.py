# The purpose of this class is to check if the directory exists and create if necessary
import os, errno


class dirCheck:

	def __init__(self):
		self.path = ""

	def checkDir(self, dir):

		self.path = dir

		if not os.path.exists(self.path):
			try:
				os.makedirs(self.path)
			except OSError as e:
				if e.errno != errno.EEXIST:
					raise	