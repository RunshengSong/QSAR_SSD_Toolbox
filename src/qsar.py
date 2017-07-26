'''
Created on Jul 26, 2017

@author: runsheng
'''

from keras.models import load_model
from sklearn.externals import joblib


class qsar():
    def __init__(self, model_path):
        self.model = self._load_model(model_path)
        
    
    def predict(self, SMILEs):
        assert self.model
        
    
    def _calculate_descriptors(self, SMILEs):
        pass

    def _load_model(self, model_path):
        return load_model(model_path)
        

if __name__ == '__main__':
    pass
