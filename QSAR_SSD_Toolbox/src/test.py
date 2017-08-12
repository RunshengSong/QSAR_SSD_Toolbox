'''
Created on Aug 5, 2017

@author: runsheng
'''
import sys
sys.path.insert(0, './src')
from scipy.stats import lognorm

from qsar import *

if __name__ == '__main__':
    SMILEs = ['CO\N=C(C(=O)OC)\c1ccccc1CO\N=C(/C)c2cccc(c2)C(F)(F)F'] # The input SMILEs must be a list

    this_q = run_all.run(SMILEs)
    
    print this_q
    
    from QSAR_SSD_Toolbox.src.ssd import ssd_generator

    this_ssd = ssd_generator()
    this_ssd.generate(this_q, 'CO\N=C(C(=O)OC)\c1ccccc1CO\N=C(/C)c2cccc(c2)C(F)(F)F', dist=lognorm, run_bootstrap=True, bootstrap_time=50, display_range=[4,10])