import numpy as np

def logSpace(vector):
		amax = np.amax(vector)
		log_vector = amax + np.log(np.exp(vector - amax))
		return log_vector


if __name__ == "__main__":

	# beta = np.array([[0.2, 0.2, 0.6],[0.1, 0.1, 0.8],[0.3, 0.4, 0.3],[0.4, 0.2, 0.4],[0.1, 0.1, 0.9],[0.2, 0.7, 0.1]])
	# qiou = np.zeros((3,3))
	V = 10000
	D = 900
	kapa = 40

	vita = np.random.rand(V,kapa) # the genes will be 180 and the topics 10.
	qiou = np.zeros((D,kapa))
	pi = np.full((1,kapa), 1.0/kapa, dtype = float)
	# coo = np.array([[1,0,1,1,1,0],[0,0,1,0,1,1],[0,1,1,1,0,1]])
	# tran_beta = beta.transpose()

	# Initialisation of coo binary matrix
	coo = np.zeros((D,V), dtype = int)
	for i in range(D):
		coo[i] = np.random.randint(2, size=V)

	# for i in range(tran_beta.shape[0]):
	# 	for j in range(coo.shape[0]):
	# 		temp = (tran_beta[i]**x[j])*(1-tran_beta[i])**(1-x[j]) 
	# 		qiou[i][j] = np.prod(temp)*0.33

	# for i in range(tran_beta.shape[0]):
	# 	qiou[i] = qiou[i]/np.sum(qiou[i])
	# # print(qiou)
	flag = False
	loop = 0
	while(flag != True):
		np.seterr(all = "warn" )
		for i in range(qiou.shape[0]): # this is actually bgc 
			# prod_vector = np.zeros((1,10), dtype = float)
			sum_vector = np.zeros((1,kapa), dtype = float)
			for j in range(qiou.shape[1]): # this is the BGC corpus
				temp = ((vita[:,j])**(coo[i]))*((1 - vita[:,j])**(1 - coo[i]))
				sum_vector[0,j] = np.log(pi[0,j]) + np.sum(np.log(logSpace(temp)))
			amax = np.amax(sum_vector)
			prod_vector = sum_vector - amax # this is the maximum value of the vector.
			reduce_factor = amax + np.log(np.sum(np.exp(prod_vector)))
					# avoid zero values.
			for j in range(qiou.shape[1]):
				if np.exp(sum_vector[0,j] - reduce_factor) < 0.0001 :
					qiou[i][j] = 0.0001
				else:
					qiou[i][j] = np.exp(sum_vector[0,j] - reduce_factor)
			qiou[i] = qiou[i]/np.sum(qiou[i])

		print("iter: {}\t qiou_min: {}".format(loop, np.amin(qiou)))
		print("iter: {}\t qiou_max: {}".format(loop, np.amax(qiou)))
		

		denominator = np.zeros((1,kapa), dtype = float)		
		numerator = np.zeros((1,kapa), dtype = float)

		for i in range(vita.shape[0]):
			for j in range(qiou.shape[1]):
				# the numerator
				num_zero = qiou[:,j]*coo[:,i]
				numerator = np.log(num_zero[ num_zero!=0])
				nmax = np.amax(numerator)
				numer_max = numerator - nmax
				num_factor = nmax + np.log(np.sum(np.exp(numer_max)))
				# the denominator
				denominator = np.log(qiou[:,j])
				dmax = np.amax(denominator)
				denom_max = denominator - dmax
				reduce_factor = dmax + np.log(np.sum(np.exp(denom_max)))
				# updating vita
				# avoid zero values
				if np.exp(num_factor - reduce_factor) <= 0.0001:
					vita[i][j] = 0.0001
				else:
					vita[i,j] = np.exp(num_factor - reduce_factor)
			# print(vita[i])
			# vita[i] = vita[i]/np.sum(vita[i])
		# update pi as well
		
		print("iter: {}\t vita_min: {}".format(loop, np.amin(vita)))
		print("iter: {}\t vita_max: {}".format(loop, np.amax(vita)))
		for i in range(qiou.shape[0]):
			pi += qiou[i]
		pi = pi/D

		if loop == 2:
			flag = True
		
		loop += 1
	print(vita[3])
	print(vita[224])