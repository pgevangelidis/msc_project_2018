# This class contains the NMF model calculations!
# Actually, it uses the scipy library and the NMF model
# However, it created the sparse matrix and computes the W and H matrices
from sklearn.decomposition import NMF
from scipy import sparse
import numpy as np

class NMF_model:
    def __init__(self, c):
        self.row_ind = []
        self.col_ind = []
        self.coordinates = c

    def setSparse(self):
        for c in self.coordinates:
            self.row_ind.append(c[0])
            self.col_ind.append(c[1])
        
        mat_coo = sparse.coo_matrix((mat_values, (row_ind, col_ind)))
        print("The sparse matrix is ready!")
        print("The number of componets are: %s\nThe number of features are: %s" % (mat_coo.shape[0], mat_coo.shape[1]))
        return mat_coo

    def modelNMF(self):
        
        sparse = self.setSparse()
        self.model = NMF(n_components=50, init='nndsvd', solver='cd')
        self.W = model.fit_transform(sparse)
        self.H = model.components_
        # Drums please....
        print('NMF is done.')

    def getMatrixW(self):
        return self.W

    def getMatrixH(self):
        return self.H
        

        