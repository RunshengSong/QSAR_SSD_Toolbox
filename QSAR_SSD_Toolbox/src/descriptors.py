'''
Created on Jul 26, 2017

@author: runsheng
'''
import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors

import warnings

class descriptor_calculator():
    def __init__(self):
        pass

    @staticmethod
    def calculate(SMILEs, filter=None):
        
        calc = Calculator(descriptors, ignore_3D=True)
        d = []
        for smi in SMILEs:
            try:
                m = Chem.MolFromSmiles(smi)
                d.append(calc(m))
            except:
                # The input SMILEs is invaild
                raise ValueError("Bad SMILEs Detected. Please Check: "+smi)
#                 warnings.warn("Bad SMILEs  Detect. Filling NA Values: "+smi)
#                 d.append(['NA'] * len(calc))
        
        d_df = pd.DataFrame(d, index=SMILEs, columns=[str(e_d) for e_d in calc.descriptors]).apply(pd.to_numeric, errors='coerce')
        if filter:
            d_df = d_df.loc[:,filter]
        
        d_df.fillna(0, inplace=True)
        return d_df.values

if __name__ == '__main__':
    # test
    pass