#This program is named "plots" but it will check if the lower bound is constantly increasing
#Furthermore the program will plot some random objects as well as the total lower bound.
import glob
import os.path
import numpy as np
import matplotlib.pyplot as plt
import re

# This class creates objects of BGC that contains the information of the lower bound progression and
# further info that will be useful to plot.
class BGCplot:
    def __init__(self):
        self.bgcName = ""
        self.lb_array = []
        self.desceding = False

    def setLBarray(self, text):
        lb_list = []
        t = []
        lb = 0.0
        for line in text:
            if 'name:' in line:
                self.bgcName = line.split()[1]
            if 'lower bound:' in line:
                if 'pre' in line:
                    temp = re.findall('\[([^[\]]*)\]', next(text))
                    t = temp.pop().split(',')
                    lb_list = list(map(float,t))
                else:
                    lb = float(next(text))
                    #print(lb)
        lb_list.remove(0.)
        lb_list.append(lb)
        self.lb_array = np.asarray(lb_list)
        # This technique converts the negative values to positive, for displaying reasons only.
        self.lb_array = (-np.amin(self.lb_array) + self.lb_array)
        self.derivativeLB()

    # this method check if the lb vector increases in the value field
    # This method is called from the setLBarrayself.
    # Inside method.
    def derivativeLB(self):
        for i in range(1,self.lb_array.shape[0]):
            if (self.lb_array[i] - self.lb_array[i-1]) < 0.0:
                self.descending = True


    def plot_LB(self):
        if (self.lb_array.size)>0:
            plt.plot(self.lb_array)
            plt.title(r'$\mathcal{\bar {LB}} - min(\mathcal{\bar {LB}})$', fontsize=18)
            plt.ylabel('Lower Bound', fontsize = 14)
            plt.xlabel('loop (#)', fontsize = 14)
            plt.legend((self.bgcName,), shadow=True, loc=(0.65, 0.05))
            plt.show()
        else:
            print('The lb array is empty.\nThe program cannot create a figure.')
##################### Class LDA ######################
# This class contains the lower bound values of all objects as a sum
######################################################
class LDAplot:
    def __init__(self):
        self.lb_array = []
        self.loop = []
        self.descending = False

    def setLBarray(self, text):
        l = []
        v = []
        for line in text:
            if not 'loop#' in line:
                t = line.split()
                l.append(int(t[0]))
                v.append(float(t[1]))
        l.reverse()
        v.reverse()
        self.loop = np.asarray(l)
        self.lb_array = np.asarray(v)
        self.lb_array = np.log(-np.amin(self.lb_array) + self.lb_array + 1)

    def derivativeLB(self):
        for i in range(1,self.lb_array.shape[0]):
            if (self.lb_array[i] - self.lb_array[i-1]) < 0.0:
                self.descending = True

    def plot_LB(self):
        if (self.lb_array.size)>0:
            plt.plot(self.lb_array)
            plt.title(r'$\log (\mathcal{\bar {LB}} - min(\mathcal{\bar {LB}}))$', fontsize=18)
            plt.ylabel('Total Lower Bound', fontsize = 14)
            plt.xlabel('loop (#)', fontsize = 14)
            plt.legend(('Total',), shadow=True, loc=(0.65, 0.05))
            plt.show()
        else:
            print('The lb array is empty.\nThe program cannot create a figure.')

##############################
#_________ Main______________#
##############################
if __name__ == '__main__':

    # The path to all BGC text files:
    bgc_path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src\bgc_objects'
    # And this path contains the lda objects
    lda_path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src\lda_object'

    # create a list with all the BGCplot objects
    bgcPlotList = []

    for filename in glob.glob(os.path.join(bgc_path, '*.txt')):
        bgc = open(filename, 'r')
        plot_obj = BGCplot()
        plot_obj.setLBarray(bgc)
        plot_obj.derivativeLB()
        if (plot_obj.desceding):
            plot_obj.plot_LB()
        bgc.close()
        bgcPlotList.append(plot_obj)

    lda_file = os.path.join(lda_path, 'totalLowerBound.txt')
    lda = open(lda_file, 'r')
    lda_plot = LDAplot()
    lda_plot.setLBarray(lda)
    lda_plot.derivativeLB()
    lda_plot.plot_LB()
    lda.close()

    print('Done')
