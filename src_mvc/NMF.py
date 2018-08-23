# This class contains the NMF model calculations!
# Actually, it uses the scipy library and the NMF model
# However, it created the sparse matrix and computes the W and H matrices
from sklearn.decomposition import NMF
from scipy import sparse
from dirCheck import *
import numpy as np

class NMF_model:
    def __init__(self):
        self.row_ind = []
        self.col_ind = []
        self.mat_values = []
        self.coordinates = []
        self.W = []
        self.H = []
        # mac path
        directory = os.path.dirname(os.path.realpath(__file__))
        self.path = directory + "/NMF_files/"
        # self.path = r'/Users/pavlos/Documents/personal/msc_project_2018/src_mvc/NMF_files/'
        # windows path
        # self.path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src_mvc\msc_project_2018\src_mvc\NMF_files'
        folder = dirCheck()
        folder.checkDir(self.path)

    def setSparse(self):
        for c in self.coordinates:
            self.row_ind.append(c[0])
            self.col_ind.append(c[1])
            self.mat_values.append(1)
        
        mat_coo = sparse.coo_matrix((self.mat_values, (self.row_ind, self.col_ind)))
        print("The sparse matrix is ready!")
        print("The number of componets are: %s\nThe number of features are: %s" % (mat_coo.shape[0], mat_coo.shape[1]))
        return mat_coo

    def modelNMF(self, c):

        self.coordinates = c
        sparse = self.setSparse()
        model = NMF(n_components=50, init='nndsvd', solver='cd')
        self.W = model.fit_transform(sparse)
        H_trans = model.components_
        self.H = np.transpose(H_trans)
        # Drums please....

        self.normaliseNMF()
        self.exportNMF()

        print('NMF is done.')
        print('W: {} H: {}'.format(self.W.shape, self.H.shape))


    def normaliseNMF(self):
        for i in range(self.H.shape[0]):
            if np.sum(self.H[i])>0:
                self.H[i] = self.H[i]/np.sum(self.H[i])

        for i in range(self.W.shape[0]):
            if (np.max(self.W[i])):
                self.W[i] = self.W[i]/np.max(self.W[i])

    def exportNMF(self):

        W_file = os.path.join(self.path, 'W_file.txt')
        H_file = os.path.join(self.path, 'H_file.txt')

        W_out = open(W_file,'w')
        H_out = open(H_file,'w')

        for i in range(self.H.shape[0]):
            H_out.write('{}'.format(self.H[i]))

        for i in range(self.W.shape[0]):
            W_out.write('{}'.format(self.W[i]))

        W_out.close()
        H_out.close()



        

        