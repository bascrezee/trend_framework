
import numpy as np
import os
import iris
import iris.pandas as ipd
import iris.plot as iplt
from scipy.stats.mstats import linregress
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf
from statsmodels.tsa.seasonal import seasonal_decompose

import copy
import matplotlib.pyplot as plt
import pprint



class TrendLims1D:
    def __init__(self,name,datatag,verbose=True,params={'alpha' : 0.05, 'trend_magnitude' : None}):
        self.name = name
        self.datatag = datatag
        self.verbose = verbose
        self.params = params
        # Initialize some empty data
        self.stats_summary = {}
        self.logbook = []
        self.__load_data__()
        self.fitted = {}
    def __load_data__(self):
        if self.datatag=='artificial':
            self.add_to_logbook("Created an artificial timeseries")
            # Create 100 years of artificial data
            timeline = pd.date_range('1/1/1901', periods=100, freq='Y')
            # Convert to time in floats, counting in decades from starting time
            timediff = timeline-pd.datetime(1901,1,1)
            decadal_time = timediff.total_seconds()/(3600.*24.*365.25*10) # Convert to units decades
            
            # Parameters for creating random test data
            mean = 1
            trend_magnitude = 0.03
            noise_magnitude = 0.05
            jump_magnitude = 0.05
            jump_start = 20
            jump_end = 40
            
            jump_axis = np.zeros_like(decadal_time)
            jump_axis[20:40] = 1
            noise_axis = np.random.uniform(-1,1,len(decadal_time))
            
            artificial_data = mean + trend_magnitude*decadal_time + jump_magnitude*jump_axis + noise_magnitude*noise_axis
            self.data_ts = pd.Series(data=artificial_data,index=timeline)
        elif os.path.splitext(self.datatag)[-1]=='.nc':
            self.add_to_logbook("Loaded datafile {0}".format(self.datatag))
            datafile = os.path.join('./test_data/',self.datatag)
            # Load the data
            self.data_cube = iris.load_cube(datafile)
            self.data_ts = iris.pandas.as_series(self.data_cube)
            self.data_ts_copy = copy.deepcopy(self.data_ts)
    def reset(self):
        self.data_ts = self.data_ts_copy
        self.logbook = self.logbook[1:] # Only preserve first line, with info on loading            
    def subset(self,timeslice):
        '''
        Subsetting over time.
        '''
        self.data_ts = self.data_ts[timeslice]
        self.add_to_logbook('Subsetted to timeperiod {0}-{1}'.format(timeslice.start,timeslice.stop))
    
    def decompose_seasonal(self,inplace=False,freq=None):
        results = seasonal_decompose(self.data_ts,model='additive',freq=freq)        
        signal = results.trend
        seasonal = results.seasonal
        residue = results.resid
        # Note that the seasonal fit, is a sum of the seasonal sign and the dataset mean
        self.fitted['seasonal'] = results.seasonal #+ self.data_ts.mean()
        if inplace:
            self.data_ts = self.data_ts-(self.fitted['seasonal'])
        else:
            return results
    
    def resample(self,target_freq):
        '''
        Resampling to the target freq
        '''
        self.data_ts = self.data_ts.resample(target_freq).mean()
        self.add_to_logbook('Resampled to {0} frequency'.format(target_freq))
    
    def get_trends(self):
        self.trend_mktest()
        self.trend_linear()

    def get_breakpoints(self):
        raise NotImplementedError
        
    def plot(self,label=None):
        self.add_to_logbook("Creating a plot with label: {0}".format(label))
        self.data_ts.plot(label=label)
        
    def trend_mktest(self):
        mk_input = self.data_ts.values
        mk_input = mk_input[~np.isnan(mk_input)]
        trend,h,p,z = self.mk_test(mk_input,alpha=self.params['alpha'])
        # Store results in stats_summary
        results = {}
        results['trend']=trend
        results['h']=h
        results['p']=p
        results['z']=z
        self.stats_summary['trend_mannkendall'] = results
        self.add_to_logbook('Calculated Mann-Kendall test')
        
    def trend_linear(self):
        xaxis = self.data_ts.index.to_julian_date()
        yaxis = self.data_ts.values
        results = linregress(xaxis,yaxis)
        self.fitted['trend_linear'] = xaxis*results.slope+results.intercept
        results_dict = {}
        results_dict['slope']=results.slope*(365.25*10) # From /day to /decade
        if results.pvalue <= self.params['alpha']:
            results_dict['trend'] = int(np.sign(results.slope))
        else:
            results_dict['trend'] = int(0)
        results_dict['pvalue']=results.pvalue
        results_dict['stderr']=results.stderr*(365.25*10)
        results_dict['slope_low'] = results_dict['slope']-results_dict['stderr']
        results_dict['slope_up'] = results_dict['slope']+results_dict['stderr']
        self.stats_summary['trend_linear'] = results_dict
        self.add_to_logbook('Calculated linear trend test')
        # We are not interested in intercept and rvalue, leave them out for now
        #results.intercept
        #results.rvalue

    # TODO
    def trend_theilsen(self):
        from scipy.stats.mstats import theilslopes

        xaxis = self.data_ts.index.to_julian_date().values
        yaxis = self.data_ts.values

        theilsen_result = theilslopes(yaxis,x=xaxis,alpha=self.params['alpha'])
        slope,intercept,slope_low,slope_up = theilsen_result
        assert(slope_low <= slope <= slope_up) # Just to be safe, check this
        slope_sign = np.sign(slope)
        if not slope_low < 0.0 < slope_up:
            trend = int(np.sign(slope))
        else:
            trend = int(0)
        
        results_dict = {}
        results_dict['trend'] = trend
        results_dict['slope'] = slope*(365.25*10) # From /day to /decade
        results_dict['slope_low'] = slope_low*(365.25*10) # From /day to /decade
        results_dict['slope_up'] = slope_up*(365.25*10) # From /day to /decade
        self.stats_summary['trend_theilsen'] = results_dict
        self.add_to_logbook('Calculated Theil-Sen slope')
        
        
    def do_trends(self):
        self.trend_mktest()
        self.trend_linear()
        self.trend_theilsen()

    def do_residual_analysis(self,fit_name=None):
        '''
        # The layout of this plot follows an example by Christoph Frei in his course 'Analysis of Weather
        # and Climate Data' at ETH Zurich
        '''
        res = self.data_ts.values-self.fitted[fit_name]

        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,squeeze=False,constrained_layout=True,figsize=(8,6))

        # ax1 QQ plot
        sm.qqplot(res/np.std(res),line='45',ax=ax1)
        ax1.set_title('Normal Q-Q Plot')
        # ax2 Histogram
        ax2.hist(res)
        ax2.set_title('Residual histogram')
        ax2.set_ylabel('# of samples')
        # testing homoscedasticity
        ax3.scatter(self.fitted[fit_name],res)
        ax3.set_title('Tukey-Anscombe')
        ax3.set_ylabel('residues')
        ax3.set_xlabel('fitted')

        #TODO this needs predicted values  y = intercept + slope*x
        nlags = 14

        autocorr = acf(res,nlags=nlags)
        lag_axis = np.array(range(nlags+1))
        
        # Now calculate ACF critical values
        # From: https://stats.stackexchange.com/questions/185425/how-to-determine-the-critical-values-of-acf
        acf_crit = np.array([1.96/np.sqrt(len(res)-x) if x!=0 else np.nan for x in lag_axis])
        ax4.plot(lag_axis,acf_crit)
        ax4.plot(lag_axis,-1*acf_crit)
        
        ax4.bar(lag_axis,autocorr)
        
        ax4.set_title('Residuals ACF')
        ax4.set_xlabel('Lag')
        ax4.set_ylabel('ACF')
        ax4.set_xticks(lag_axis)
        ax4.set_xticklabels(lag_axis)
        return f

    def mk_test(self,x, alpha=0.05):
        """
        This function is derived from code originally posted by Sat Kumar Tomer
        (satkumartomer@gmail.com)
        See also: http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm
        The purpose of the Mann-Kendall (MK) test (Mann 1945, Kendall 1975, Gilbert
        1987) is to statistically assess if there is a monotonic upward or downward
        trend of the variable of interest over time. A monotonic upward (downward)
        trend means that the variable consistently increases (decreases) through
        time, but the trend may or may not be linear. The MK test can be used in
        place of a parametric linear regression analysis, which can be used to test
        if the slope of the estimated linear regression line is different from
        zero. The regression analysis requires that the residuals from the fitted
        regression line be normally distributed; an assumption not required by the
        MK test, that is, the MK test is a non-parametric (distribution-free) test.
        Hirsch, Slack and Smith (1982, page 107) indicate that the MK test is best
        viewed as an exploratory analysis and is most appropriately used to
        identify stations where changes are significant or of large magnitude and
        to quantify these findings.
        
        #############################################################################
        MIT License
        Copyright (c) 2017 Michael Schramm
        
        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:
        
        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.
        
        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.
        #############################################################################    

        
        Input:
            x:   a vector of data
            alpha: significance level (0.05 default)
        Output:
            trend: tells the trend (increasing, decreasing or no trend)
            h: True (if trend is present) or False (if trend is absence)
            p: p value of the significance test
            z: normalized test statistics
        Examples
        --------
          >>> x = np.random.rand(100)
          >>> trend,h,p,z = mk_test(x,0.05)
        """
        n = len(x)

        # calculate S
        s = 0
        for k in range(n-1):
            for j in range(k+1, n):
                s += np.sign(x[j] - x[k])

        # calculate the unique data
        unique_x = np.unique(x)
        g = len(unique_x)

        # calculate the var(s)
        if n == g:  # there is no tie
            var_s = (n*(n-1)*(2*n+5))/18
        else:  # there are some ties in data
            tp = np.zeros(unique_x.shape)
            for i in range(len(unique_x)):
                tp[i] = sum(x == unique_x[i])
            var_s = (n*(n-1)*(2*n+5) - np.sum(tp*(tp-1)*(2*tp+5)))/18

        if s > 0:
            z = (s - 1)/np.sqrt(var_s)
        elif s < 0:
            z = (s + 1)/np.sqrt(var_s)
        else: # s == 0:
            z = 0

        # calculate the p_value
        p = 2*(1-norm.cdf(abs(z)))  # two tail test
        h = abs(z) > norm.ppf(1-alpha/2)

        if (z < 0) and h:
            trend = -1
        elif (z > 0) and h:
            trend = 1
        else:
            trend = 0
        return trend, h, p, z
    def add_to_logbook(self,logmessage):
        if self.verbose:
            print(logmessage)
        self.logbook.append(logmessage+'\n')
    def print_logbook(self):
        print('\n'.join(self.logbook))
    def print_stats(self):
        pprint.pprint(self.stats_summary)
        
def assess_trend_consistency(trendobject):
    trend_consistency_score = 0
    all_trends = [trendobject.stats_summary[trendname]['trend'] for trendname in ['trend_linear','trend_mannkendall','trend_theilsen']]
    trend_text = {
        -1 : 'a negative trend',
        0  : 'no trend',
        1 :  'a positive trend'
    }
    if len(set(all_trends))==1:
        print("The different trend tests are fully consistent and indicate {0}.".format(trend_text[list(set(all_trends))[0]]))
        trend_consistency_score += 2
    elif len(set(all_trends))==3:
        print("The different trend tests are fully in consistent")
    elif len(set(all_trends))==2:
        most_common = max(set(all_trends), key=all_trends.count)
        print("Two out of three trend tests indicate {0}.".format(trend_text[most_common]))
        trend_consistency_score += 1 # Dummy statement, just for completion of the scoring
    slope_consistency_score = 0
    name_a,name_b='trend_linear','trend_theilsen'
    slope_a_within_range_b = trendobject.stats_summary[name_b]['slope_low'] <= trendobject.stats_summary[name_a]['slope'] <= trendobject.stats_summary[name_b]['slope_up']
    slope_b_within_range_a = trendobject.stats_summary[name_a]['slope_low'] <= trendobject.stats_summary[name_b]['slope'] <= trendobject.stats_summary[name_a]['slope_up']
    overlap_of_range = max(trendobject.stats_summary[name_a]['slope_low'],trendobject.stats_summary[name_a]['slope_up']) <= min(trendobject.stats_summary[name_b]['slope_low'],trendobject.stats_summary[name_b]['slope_up'])
    if slope_a_within_range_b and slope_b_within_range_a:
        print("The two different slopes are fully consistent")
        slope_consistency_score += 3
    elif slope_a_within_range_b or slope_b_within_range_a:
        print("The two different slopes are partially consistent")
        slope_consistency_score += 2
    else:
        if overlap_of_range:
            print("The two different slope ranges marginally overlap.")
            slope_consistency_score += 1
        else:
            print("The two different slopes are very inconsistent")
            slope_consistency_score += 0 # Dummy statement, just for completion of the scoring
    trend_score = trend_consistency_score+slope_consistency_score
    print("Final trend score: {0}".format(trend_score))
    return trend_score
        
