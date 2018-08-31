class lb_check:
    
	def __init__(self):
		self.partA = []
		self.partBE = []
		self.partC = []
		self.partD = []

	def check_curve(self, factor, name):
		# the parameter factor is a list of float numbers.
		if not factor:
			return "The list is empty"
		else:
			for i in range (0, len(factor)-1 ):
				pre = factor[i]
				post = factor[i+1]
				if (post - pre) < 0:
					print("In {} iteration the curve is decreasing in part {}\nPre: {} Post: {}".format(i, name, pre, post))
			print("The process is done")
				





