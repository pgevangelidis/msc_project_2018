# this is a test script to check if my objects can be modified in other classes

import numpy as np
import scipy.special
import scipy.integrate as integrate
import glob
import os.path

if __name__ == '__main__':

	fooList = []
	beta = np.array(([0.2, 0.2, 0.6],[0.1, 0.1, 0.8],[0.3, 0.4, 0.3],[0.4, 0.2, 0.4]))
	w = np.array(([1,1,0,1],[0,0,1,0],[0,1,1,1]))
	qiou = np.zeros((3,3), dtype = float)
	p = np.full((3,1), 1/3)
	betaT = beta.transpose()

	print(betaT)
	print(betaT.shape)
	print(w.shape)
	for j in range(w.shape[0]):
		for i in range(betaT.shape[0]):
			temp = p[i]*(betaT[i]**w[j])*(1 - betaT[i])**(1 - w[j])
			qiou[i][j] = np.prod(temp)

	for i in range(qiou.shape[0]):
		qiou[i] = qiou[i]/np.sum(qiou[i])
		print(qiou[i])

	num = np.zeros((1,3))
	for i in range(w.shape[0]):
		num += qiou[i]
	print(num)

	wT = w.transpose()
	print(wT)
	# qiouT = qiou.transpose()
	# print(qiouT)

	beta = np.prod(wT, qiou)
	print(beta)

	print("Done")
	


	### now update the beta!
