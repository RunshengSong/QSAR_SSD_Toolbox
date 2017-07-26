'''
Created on Jul 26, 2017

@author: runsheng
'''
import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors

class descriptor_calculator():
    def __init__(self, filter=None):
        self.filter = filter
        self.calc = Calculator(descriptors, ignore_3D=True)
    
    def calculate(self, SMILEs):
        d = []
        for smi in SMILEs:
            try:
                m = Chem.MolFromSmiles(smi)
                d.append(self.calc(m))
            except:
                d.append(['NA'] * len(self.calc))
        
        self.d_df = pd.DataFrame(d, index= SMILEs, columns=[str(e_d) for e_d in self.calc.descriptors()]).apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all', inplace=True)
        
        if self.filter:
            self.d_df = self.d_df[self.filter]
        
        return self.d_df

if __name__ == '__main__':
    # test
    pass