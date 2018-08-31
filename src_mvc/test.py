import numpy as np

if __name__ == "__main__":

	# beta = np.array([[0.2, 0.2, 0.6],[0.1, 0.1, 0.8],[0.3, 0.4, 0.3],[0.4, 0.2, 0.4],[0.1, 0.1, 0.9],[0.2, 0.7, 0.1]])
	# qiou = np.zeros((3,3))
	V = 1800
	D = 75
	kapa = 10

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

	for i in range(qiou.shape[0]): # this is actually bgc 
		# prod_vector = np.zeros((1,10), dtype = float)
		sum_vector = np.zeros((1,10), dtype = float)
		for j in range(qiou.shape[1]): # this is the BGC corpus
			temp = ((vita[:,j])**(coo[i]))*((1 - vita[:,j])**(1 - coo[i]))
			# print("{}\tpi {}\tmultiple pi {}\tlog {}".format(np.prod(temp), pi[0,j], pi[0,j]*np.prod(temp), np.log(np.prod(temp))))
			# prod_vector[0,j] = pi[0,j]*np.prod(temp) # the dimensions are D x k. 
			sum_vector[0,j] = np.log(pi[0,j]) + np.sum(np.log(temp))
		# print("Pre update qiou: {}".format(qiou[i]))
		# print(sum_vector)
		reduce_factor = np.sum(np.exp(sum_vector))
		print(sum_vector)
		print(reduce_factor)
		for j in range(qiou.shape[1]):
			qiou[i][j] = np.exp(sum_vector[0,j])/reduce_factor
		# print(qiou[i])
	print(np.amin(qiou))
	print(np.amax(qiou))

