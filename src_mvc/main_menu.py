# This is the main file
# the main is responsible to read all the .gbk files
# create BGC objects and store all the necessary info
# create two dictionaries for the sparse matrix
# And initialise the LDA


#### packages/libraries/ etc
from Bio import *
from BGC import *
from BGC_Dictionary import *
from LDAmodel import *
from storeBGC import *
from showMenu import *
from Controller import *
import glob
import os.path
import time
import copy

bgcList = []

if __name__ == '__main__':

	print('-----------------------------')
	print('---- smCOG Decomposition ----')
	print('-----------------------------')
    # Those are the models
	nmf_flag = False
	lda_original_flag = False
	lda_binary_flag = False
    # This flag directs the program to the pickle class
	load_bgc_pickle = False
	# If the user decides to exit the program the flag will turn to True
	exitCondition = False
	# this is the reference for the showMenu class
	menuToString = showMenu()
	# this is the reference for the controller
	controllerMenu = Controller()

	while controllerMenu.exitCond != True :
		
		# show the menu to the user by using the menu() method
		menuToString.menu()


		flag = False
		while (flag!=True) :
			try:
				mode = int(raw_input('Enter your input:'))
				if(mode>0 or mode<7):
					flag = True
			except ValueError:
				print("Not a number")

		# create the controller who will control the mode value
		controllerMenu.defineMode(mode)

	print('\nThank you. Good Bye.\n')

