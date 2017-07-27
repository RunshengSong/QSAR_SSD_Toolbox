'''
Created on Jul 26, 2017

Pipeline to create Species Sensitivity Distribution Curves

@author: runsheng
'''

import pandas as pd
import numpy as np
from scipy.stats import lognorm
from scipy.stats import norm
import matplotlib.pyplot as plt



class ssd_generator():
    def __init__(self):
        pass
    
    def generate(self, 
                 ssd_df, 
                 dist=lognorm, 
                 run_bootstrap=True,
                 bootstrap_time=1000,
                 display=True):
        '''
        Input: pandas dataframe of Species_Name: LC50 Values
        Output: SSD Curves
        '''
        ssd_df, frac = self._frac(ssd_df, fraction=0.5)
        shape, loc, scale = dist.fit(ssd_df['val'], floc=0)
        
        newx = pow(np.linspace(np.log10(0.01), np.log10(max(ssd_df['val'])), 1000), 10)
        
        average_fit = dist.cdf(newx, s=shape, loc=0, scale=scale)
        
        if display:
            # plot
            fig, ax = plt.subplots()
            if run_bootstrap:
                this_boots = bootstrap()
                boots_res = this_boots.run_boots(ssd_df, newx, shape, scale, dist=dist, times=bootstrap_time)
                upper, lower = this_boots.get_upper_lower(boots_res, upper=97.5, lower=2.5)
                for i in range(len(boots_res)):
                    ax.plot(newx, boots_res[i], c='steelblue', alpha=0.2, zorder=0)
                # plot the upper and lower CI
                ax.plot(newx, upper,':', c='black')
                ax.plot(newx, lower,':', c='black')
            
            #plot the average curve
            ax.plot(newx, average_fit,color='red')
            
            #plot the points
            ax.scatter(ssd_df['val'], frac, c='black', zorder=10)
            
            # annotation
            for i, txt in enumerate(ssd_df.index):
                ax.annotate(txt, (ssd_df['val'].values[i], frac[i]))
 
            # configeration
            plt.xscale('log')
            plt.xlim([min(ssd_df['val']) - min(ssd_df['val']*0.1), max(ssd_df['val']) + max(ssd_df['val'])*0.1])
            plt.grid()
            
            plt.show()
        
        return shape, loc, scale
            
    
    def _frac(self, df, fraction=0.5):
        df = df.sort(['val'])
        frac = self._ppoints(df['val'], fraction)
        return df, frac
        
    def _ppoints(self, n, a):
        """ numpy analogue or `R`'s `ppoints` function
        see details at http://stat.ethz.ch/R-manual/R-patched/library/stats/html/ppoints.html 
        :param n: array type or number
        """
        try:
            n = np.float(len(n))
        except TypeError:
            n = np.float(n)
        return (np.arange(n) + 1 - a)/(n + 1 - 2*a)
    

class bootstrap():
    def __init__(self):
        pass
    
    def run_boots(self, df, newx, shape, scale, dist=lognorm, times=1000):
        res = []
        for i in range(times):
            this_boot = self._boots(df, newx, shape, scale, dist=dist)
            res.append(this_boot)
        return np.asarray(res)
    
    def get_upper_lower(self, boots_res, upper=97.5, lower=2.5):
        return (np.percentile(boots_res, upper, axis=0), (np.percentile(boots_res, lower, axis=0)))
        
    def _boots(self, df, newx, shape, scale, dist=lognorm):
        xr = lognorm.rvs(size=len(df['val']), s=shape, loc=0, scale=scale)
        this_shape, this_loc, this_scale = lognorm.fit(xr, floc=0)
        this_fit = dist.cdf(newx, s=this_shape, loc=0, scale=this_scale)
        
        return list(this_fit)
    
    
        
if __name__ == '__main__':
    # test
    import io
    import requests
    
    url="https://raw.githubusercontent.com/EDiLD/r-ed/master/post_ssd/ssd_data.csv"
    s=requests.get(url).content
    c=pd.read_csv(io.StringIO(s.decode('utf-8')))
    c.set_index('species',inplace=True)
    
    this_generator = ssd_generator()
    this_generator.generate(c, dist=lognorm, run_bootstrap=True, bootstrap_time=1000, display=True)
    
    