'''
Created on Jul 26, 2017

@author: runsheng
'''
import os
import json
import numpy as np
import pandas as pd
from scipy.spatial import distance
from collections import defaultdict

from keras.models import load_model
from sklearn.externals import joblib

from descriptors import descriptor_calculator

class qsar():
    def __init__(self, model_name):
        self.model = self._load_model(model_name)
        self.scalar = self._load_scalar(model_name)
        self.filter = self._load_filter(model_name)
        self.ad_dict = self._load_ad(model_name)
     
    def predict(self, SMILEs):
        SMILEs = [SMILEs] # temp solution...
        
        if len(SMILEs) > 1:
            raise ValueError("Only accept sinlge SMILEs") 
        
        descriptors = self.scalar.transform(descriptor_calculator.calculate(SMILEs, self.filter))
        inside_ad, error_bars = self._calculate_ad(descriptors, self.ad_dict)
        prediction = self.model.predict(descriptors)
        prediction_higher, prediction_lower = self._get_range(prediction[0], error_bars[0])
        
        return prediction[0][0], inside_ad[0], error_bars[0], prediction_higher[0], prediction_lower[0]
    
    def _calculate_ad(self, descriptors, ad_dict):
        cut_off_threshold = ad_dict['cut off']
        error_inside = ad_dict['error inside']
        error_outside = ad_dict['error outside']
        c = ad_dict['centroid']
        
        distance_to_c = [distance.euclidean(i, c) for i in descriptors]
        inside_ad = [1 if d < cut_off_threshold else 0 for d in distance_to_c]
        error_bars = [error_inside if k == 1 else error_outside for k in inside_ad]
        
        return inside_ad, error_bars
    
    def _get_range(self, prediction, error_bars):
        return prediction + prediction*error_bars, prediction - prediction*error_bars
    
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
    
    def _load_ad(self, model_name):
        cur_path = os.path.dirname(__file__)
        with open(cur_path+'/../models/'+model_name+'/ad.json') as ad_data:
            ad_dict = json.load(ad_data)
        
        return ad_dict
    
class run_all():
    @staticmethod
    def run(SMILEs):
        cur_path = os.path.dirname(__file__)
        all_models = [d for d in os.listdir(cur_path+'/../models') if os.path.isdir(os.path.join(cur_path+'/../models', d))]
        species = []
        all_p = defaultdict(list)
 
        for each_model in all_models:
            print each_model
            species.append(each_model)
            this_qsar = qsar(each_model)
            this_p, this_inside, this_error, this_higher, this_lower = this_qsar.predict(SMILEs)
            all_p[each_model] = [this_p, this_inside, this_error, this_higher, this_lower]

        df = pd.DataFrame.from_dict(all_p, orient='index')
        df.columns = ['Prediction', 'Inside AD', 'Prediction Error', 'Prediction Upper', 'Prediction Lower']
        return df
    
if __name__ == '__main__':
    pass
#     # test
# #     model_name = "Oncorhynchus Mykiss"
# #        
# #     q = qsar(model_name=model_name)
# #     print q.predict(['CCC', 'C[N+](C)(C)CCCl.[Cl-]','ClC(Cl)(Cl)C(Cl)(Cl)Cl','c1ccccc1','ClC(Cl)(Cl)SN1C(=O)c2ccccc2C1=O','CCSC(=O)N(CC(C)C)CC(C)C'])
#     print run_all.run(['CCC'])
   
    