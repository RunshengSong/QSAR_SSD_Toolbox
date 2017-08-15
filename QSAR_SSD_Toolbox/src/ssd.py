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
                 run_error =True,
                 run_bootstrap=True,
                 bootstrap_time=1000,
                 display_range=[0,100]):
        '''
        Input: pandas dataframe of Species_Name: LC50 Values
        Output: SSD Curves
        '''
        df_mean, frac_mean, df_high, frac_high, df_low, frac_low = self._frac(ssd_df, fraction=0.5) # dataframe is sorted here
        
        shape_mean, loc_mean, scale_mean = dist.fit(df_mean['Prediction'], floc=0)
        shape_lower, loc_lower, scale_lower = dist.fit(df_high['Prediction Lower'], floc=0)
        shape_upper, loc_upper, scale_upper = dist.fit(df_low['Prediction Upper'], floc=0)
        
        newx = pow(np.linspace(np.log10(0.01), np.log10(max(df_mean['Prediction'])), 1000), 10)
        
        fit_mean = dist.cdf(newx, s=shape_mean, loc=0, scale=scale_mean)
        fit_lower = dist.cdf(newx, s=shape_lower, loc=0, scale=scale_lower)
        fit_upper = dist.cdf(newx, s=shape_upper, loc=0, scale=scale_upper)
        
        # plot
        fig, ax = plt.subplots(figsize=(11,10))
        if run_bootstrap:
            this_boots = bootstrap()
            boots_res = this_boots.run_boots(df_mean, newx, shape_mean, scale_mean, dist=dist, times=bootstrap_time)
            upper, lower = this_boots.get_upper_lower(boots_res, upper=97.5, lower=2.5)
            for i in range(len(boots_res)):
                ax.plot(newx, boots_res[i], c='steelblue', alpha=0.2, zorder=0)
            # plot the upper and lower CI
            ax.plot(newx, upper,':', c='black', label='Upper and Lower Bound of Bootstrapping SSD')
            ax.plot(newx, lower,':', c='black')
            
        #plot the average curve
        ax.plot(newx, fit_mean,color='red',label='Average SSD')
        if run_error:
            ax.plot(newx, fit_lower,color='lightsalmon',label='Lower Bound')
            ax.plot(newx, fit_upper,color='lightsalmon',label='Upper Bound')
             
        #plot the points
        ax.scatter(df_mean['Prediction'], frac_mean, c='black', zorder=10)
        if run_error:
            ax.scatter(df_low['Prediction Lower'], frac_low, marker='x',c='gray', zorder=10)
            ax.scatter(df_high['Prediction Upper'], frac_high, marker='x', c='gray', zorder=10)
            
            
        # annotation
        for i, txt in enumerate(df_mean.index):
            ax.annotate(txt, (df_mean['Prediction'].values[i], frac_mean[i]))
        if run_error:
            for i, txt in enumerate(df_low.index):
                ax.annotate(txt, (df_low['Prediction Lower'].values[i]-0.15, frac_low[i]), size='small',color='gray')
                
            for i, txt in enumerate(df_high.index):    
                ax.annotate(txt, (df_high['Prediction Upper'].values[i], frac_high[i]-0.03), size='small',color='gray')
            
        # configeration
        plt.xscale('log')
        plt.xlim(display_range)
        plt.grid()
        plt.xlabel('Concentration (mol/L)')
        plt.ylabel('Cumulative Probability')
        plt.legend(loc='upper left')
            
        plt.show()
        return df_mean, df_high, df_low, fig
            
    def _frac(self, df, fraction=0.5):
        df_mean = df.sort_values(['Prediction'])
        frac_mean = self._ppoints(df['Prediction'], fraction)
        
        df_high = df.sort_values(['Prediction Upper'])
        frac_high = self._ppoints(df['Prediction Upper'], fraction)
        
        df_low = df.sort_values(['Prediction Lower'])
        frac_low = self._ppoints(df['Prediction Lower'], fraction)
        
        return df_mean, frac_mean, df_high, frac_high, df_low, frac_low
       
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
    
    def run_boots(self, df,newx, shape, scale, dist=lognorm, times=1000):
        res = []
        for i in range(times):
            this_boot = self._boots(df,newx, shape, scale, dist=dist)
            res.append(this_boot)
        return np.asarray(res)
    
    def get_upper_lower(self, boots_res, upper=97.5, lower=2.5):
        return (np.percentile(boots_res, upper, axis=0), (np.percentile(boots_res, lower, axis=0)))
        
    def _boots(self, df, newx, shape, scale, dist=lognorm):
        xr = lognorm.rvs(size=len(df['Prediction']), s=shape, loc=0, scale=scale)
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
    
    