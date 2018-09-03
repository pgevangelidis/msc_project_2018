import glob
import os.path
import time
import copy
import re
import numpy as np
import matplotlib.pyplot as plt
from dirCheck import *

if __name__ == "__main__":

	read_path = r"/Users/pavlos/Documents/personal/msc_project_2018/src_mvc/lda_objects_binary"
	lb_binary_file = "totalLowerBound.txt"
	a_path = os.path.join(read_path, lb_binary_file)
	file_handler = open(a_path, 'r')

	lb_list = []
	next(file_handler)
	for line in file_handler:
		lb_list.append(float(line.split()[1]))
	lb_list.reverse()

	lb_array = np.asarray(lb_list)
	# lb_array = np.log(-np.amin(lb_array) + lb_array + 0.0001)
	print(np.log(lb_array - np.amin(lb_array) + 1.001))
	# lb_array = np.amax(lb_array) + np.log(np.exp(lb_array - np.amax(lb_array)))
	lb_array = np.log(lb_array - np.amin(lb_array) + 1.0001)

	if (len(lb_array))>0:
		plt.plot(lb_array)
		plt.title(r'$\log (\mathcal{\bar {LB}} - min(\mathcal{\bar {LB}}))$', fontsize=18)
		plt.ylabel('Total Lower Bound', fontsize = 14)
		plt.xlabel('loop (#)', fontsize = 14)
		plt.legend(('Total',), shadow=True, loc=(0.65, 0.05))
		plt.show()
	else:
		print('The lb array is empty.\nThe program cannot create a figure.')