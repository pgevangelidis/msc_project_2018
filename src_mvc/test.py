import numpy as np

if __name__ == "__main__":

	beta = np.array([[0.2, 0.2, 0.6],[0.1, 0.1, 0.8],[0.3, 0.4, 0.3],[0.4, 0.2, 0.4]])
	qiou = np.zeros((3,3))
	
	x = np.array([[1,0,1,1],[0,0,1,0],[0,1,1,1]])
	tran_beta = beta.transpose()

	print(beta)
	print(beta[0])
	print(beta[:,0])

	for i in range(tran_beta.shape[0]):
		for j in range(x.shape[0]):
			temp = (tran_beta[i]**x[j])*(1-tran_beta[i])**(1-x[j]) 
			qiou[i][j] = np.prod(temp)*0.33
			# print("\niteration {}\n".format(i))
			# print("beta transpose: {}\n x vector: {}".format(tran_beta[i], x[j]))
			# print((tran_beta[i]**x[j]))
			# print((1-tran_beta[i])**(1-x[j]))
			# print("the actual vector is:\n{}\n*******".format(temp))

	for i in range(tran_beta.shape[0]):
		qiou[i] = qiou[i]/np.sum(qiou[i])
	# print(qiou)

	print(np.log(0.0000001))
	# print(beta.transpose()[0])

