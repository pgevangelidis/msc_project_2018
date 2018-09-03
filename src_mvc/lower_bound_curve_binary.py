import glob
import os.path
import time
import copy
import re
import numpy as np
from dirCheck import *

if __name__ == "__main__":

	read_path = r"C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\lda_objects_binary"
	lb_binary_file = "totalLowerBound"
	a_path = os.path.join(read_path, lb_binary_file)
	file_handler = open(a_path, 'r')

	lb_array = []

	if (len(lb_array))>0:
            plt.plot(lb_array)
            plt.title(r'$\log (\mathcal{\bar {LB}} - min(\mathcal{\bar {LB}}))$', fontsize=18)
            plt.ylabel('Total Lower Bound', fontsize = 14)
            plt.xlabel('loop (#)', fontsize = 14)
            plt.legend(('Total',), shadow=True, loc=(0.65, 0.05))
            plt.show()
        else:
            print('The lb array is empty.\nThe program cannot create a figure.')