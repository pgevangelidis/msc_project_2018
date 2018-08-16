# Class foo_doo
from foo import *

class foo_doo:
	def __init__(self, g):
		self.parameter = g

	def testCounter(self, x):
		x.counter += 5

	def testName(self, x, s):
		x.name = s