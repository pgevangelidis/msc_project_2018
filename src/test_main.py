# this is a test script to check if my objects can be modified in other classes
from foo import *
from foo_doo import *
import numpy as np
import scipy.special
import scipy.integrate as integrate
import glob
import os.path

if __name__ == '__main__':

	# fooList = []
	# for i in range(6):
	#
	# 	x = foo(1)
	# 	if i==2 or i==4:
	# 		x.setCounter()
	# 	fooList.append(x)
	# 	print("the name: {} and the counter: {}".format(x.name, x.counter))
	# 	print(x)
	#
	# z = foo_doo(2)
	# for y in fooList:
	# 	z.testCounter(y)
	# 	z.testName(y, "test name")
	# 	print("the NEW name: {} & counter: {}".format(y.name, y.counter))
	# 	print(y)
	#
	# ### Now I will modify the x object in
	# print('I small test:\n {}'.format(fooList[1].name))

	# gamma = np.full((1,50), 1.02)
	# phi = np.full((1,50), 0.02)
	# vita = np.random.rand(1, 50)
	# print('***Initial vectors***\n')
	# print('phi:\n{}\ngamma:\n{}\nvita:\n{}\n'.format(phi, gamma, vita))
	# for i in range(5):
	# 	factor_phi = scipy.special.digamma(np.sum(gamma))
	# 	for j in range(50):
	# 		factor = scipy.special.digamma(gamma[0][j])
	# 		phi[0][j] = vita[0][j]*np.exp(factor - factor_phi)
	# 	phi = np.exp(phi - np.amax(phi))
	# 	gamma = 1 + phi
	# 	print('***\nLoop number: {}\n'.format(i))
	# 	print('gamma {}\nsum gamma: {}\n'.format(gamma, np.sum(gamma)))
	# 	print('phi: {}\n'.format(phi))
	# 	print('factor_phi: {}\tfactor: {}\nnew phi:{}\n'.format(factor_phi, factor, phi))

	bgc_path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src\bgc_objects_test'
	for filename in glob.glob(os.path.join(bgc_path, '*.txt')):
		bgcName = ""
		noGenes = 0
		geneList = []
		gam = []
		bgc = open(filename, 'r')
		for line in bgc:
			if 'name:' in line:
				bgcName = line.split()[1]
				print(bgcName)
			# if 'genes:' in line:
			# 	noGenes = int(next(bgc))
			# 	for i in range(noGenes):
			# 		t = next(bgc).split()
			# 		if int(t[2]) > 1:
			# 			print(bgcName)
			# 			break
			# 	break
			if 'gamma:' in line:
				t = next(bgc).split()
				t = list(map(float, t))
				gam = np.asarray(t)
				# gam = gam/np.sum(gam)
				# gam = gam/np.linalg.norm(gam) # This is a standard normalise method but it can not satisfy the sum to 1
				# gam = gam/(np.sum(gam))


	# partD = (-1)*(np.log(scipy.special.gamma(np.sum(gam))))
	# print(np.log(scipy.special.gamma(np.sum(gam)-100.0)))
	# print('loop 0 part D{}'.format(partD))
	# for i in range(gam.shape[0]):
	# 	partD += np.log(scipy.special.gamma(gam[i]))
	# 	partD -= (gam[i] - 1)*(scipy.special.digamma(gam[i]) - scipy.special.digamma(np.sum(gam)))
	# 	# print('loop {} part D{}'.format(i,partD))

	# for i in range(150,200,10):
	# 	print('sum {} digamma {}'.format(i, np.log(scipy.special.polygamma(0,i))))
	#
	vita_test = np.ones((5,10), dtype=float)
	half_vita = np.full((3,10), 0.759)
	quart_vita = np.full((1,50), 1.741)

	big_vita = np.vstack((vita_test, half_vita))
	# big_vita = np.vstack((big_vita, quart_vita))


	print('gamma 1st half: {}\ngamma 2nd half: {}'.format(scipy.special.gammaln(np.sum(quart_vita[0][0:9])), scipy.special.gammaln(np.sum(quart_vita[0][10:19]))))
	print('gamma 1st: {}\ngamma 2nd: {}'.format(quart_vita[0][0:4], quart_vita[0][5:9]))
	print('gamma: {}'.format(scipy.special.gammaln(np.sum(quart_vita))))

	der = scipy.special.polygamma(0, np.sum(quart_vita))
	# print(gam)
	print("Done")
