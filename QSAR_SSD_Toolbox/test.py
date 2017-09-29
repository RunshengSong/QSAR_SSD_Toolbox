'''
Created on Aug 5, 2017

@author: runsheng
'''
import sys

from scipy.stats import lognorm

from src.qsar import *

if __name__ == '__main__':
    SMILEs = ['C=CC=O'] # The input SMILEs must be a list

#     this_q = qsar("Americamysis bahia")
#     print this_q.predict(SMILEs)
    this_q = run_all.run_all_species(SMILEs)
#     
#     print this_q.predict(SMILEs)
    
    from QSAR_SSD_Toolbox.src.ssd import ssd_generator

    this_ssd = ssd_generator()
    print this_q
    raw_input()
    df_mean, df_high, df_low, this_fig = this_ssd.generate(this_q, dist=lognorm, run_error=True,run_bootstrap=False, bootstrap_time=500, display_range=[1, 10000])
# 
#     print df_mean
#     print df_high
#     print df_low
    