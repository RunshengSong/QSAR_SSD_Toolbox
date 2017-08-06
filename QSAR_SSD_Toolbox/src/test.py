'''
Created on Aug 5, 2017

@author: runsheng
'''
import sys
sys.path.insert(0, './src')
from scipy.stats import lognorm

from qsar import *

if __name__ == '__main__':
    SMILEs = ['CCCC'] # The input SMILEs must be a list

    this_q = run_all.run(SMILEs)
    
    from QSAR_SSD_Toolbox.src.ssd import ssd_generator

    this_ssd = ssd_generator()
    this_ssd.generate(this_q, dist=lognorm, run_bootstrap=True, bootstrap_time=100, display=True)