'''
Created on Jul 26, 2017

@author: runsheng
'''
import os

from keras.models import load_model
from sklearn.externals import joblib

from descriptors import descriptor_calculator


class qsar():
    def __init__(self, model_name):
        self.model = self._load_model(model_name)
        self.scalar = self._load_scalar(model_name)
        self.filter = self._load_filter(model_name)
     
    def predict(self, SMILEs):
        assert self.model
        descriptors = self.scalar.transform(descriptor_calculator.calculate(SMILEs, self.filter))
        
        prediction = self.model.predict(descriptors)
        return prediction

    def _load_model(self, model_name):
        cur_path = os.path.dirname(__file__)
        return load_model(cur_path+'/../models/'+model_name+'/model.h5')
    
    def _load_scalar(self,model_name):
        cur_path = os.path.dirname(__file__)
        return joblib.load(cur_path+'/../models/'+model_name+'/scaler.pkl')
    
    def _load_filter(self,model_name):
        cur_path = os.path.dirname(__file__)
        with open(cur_path+'/../models/'+model_name+'/filter.txt') as myfile:
            content = myfile.readlines()
        return [x.strip() for x in content]
        

if __name__ == '__main__':
    pass
    # test
#     model_name = "Oncorhynchus_Mykiss"
#       
#     q = qsar(model_name=model_name)
#     print q.predict(['CCC', 'C[N+](C)(C)CCCl.[Cl-]','ClC(Cl)(Cl)C(Cl)(Cl)Cl','c1ccccc1','ClC(Cl)(Cl)SN1C(=O)c2ccccc2C1=O','CCSC(=O)N(CC(C)C)CC(C)C'])
#     
    