#This program is named "plots" but it will check if the lower bound is constantly increasing
#Furthermore the program will plot some random objects as well as the total lower bound.
import glob
import os.path
import numpy as np
import matplotlib.pyplot as plt
import re
import copy

# This class creates objects of BGC that contains the information of the lower bound progression and
# further info that will be useful to plot.
class BGCplot:
    
    def __init__(self):
        self.bgcName = ""
        self.lb_array = []
        self.desceding = False
        self.loop = []
        self.bgcList = []

    def createPlot(self, lda, bgcList):
        # first I will iterate the total lower bound pre except the zero point
        self.bgcList = copy.deepcopy(bgcList)
        lb_list = []

        for i in range(1, len(lda.totalLBound_pre)):
            lb_list.append(lda.totalLBound_pre[i])
        lb_list.append(lda.totalLBound)

        self.lb_array = np.asarray(lb_list)
        # This technique converts the negative values to positive, for displaying reasons only.
        self.lb_array = np.log(-np.amin(self.lb_array) + self.lb_array)

        self.derivativeLB()
        self.plot_LB()
        # this method check if the lb vector increases in the value field
        # This method is called from the setLBarrayself.
        # Inside method.
    def derivativeLB(self):
        for i in range(1,self.lb_array.shape[0]):
            if (self.lb_array[i] - self.lb_array[i-1]) < 0.0:
                self.descending = True

    def plot_LB(self):
        if (len(self.lb_array))>0:
            plt.plot(self.lb_array)
            plt.title(r'$\log (\mathcal{\bar {LB}} - min(\mathcal{\bar {LB}}))$', fontsize=18)
            plt.ylabel('Total Lower Bound', fontsize = 14)
            plt.xlabel('loop (#)', fontsize = 14)
            plt.legend(('Total',), shadow=True, loc=(0.65, 0.05))
            plt.show()
        else:
            print('The lb array is empty.\nThe program cannot create a figure.')
