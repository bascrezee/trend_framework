import pandas as pd
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

# rpy2 imports
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from collections import OrderedDict

import copy
import matplotlib.pyplot as plt
import pprint

# A dictionary for converting long names to short names
breakpoint_test_names = {
    'Standard Normal Homogeneity Test (SNHT)' : 'SNHT',
    "Pettitt's test for single change-point detection" : 'Pettitt',
    'Buishand range test' : 'Buishand R',
    'Buishand U test' : 'Buishand U',
}

def rvector_to_pydict(rvector):
    '''Conversion of R vector (a Python representation of an R object, module rpy2) to a Python dictionary
    
    This solution has been taken from: https://stackoverflow.com/a/24153569/2164368
    
    Parameters
    -----------
    rvector :  rvector 
        The rvector object to be converted
        
    Returns
    --------
    dict
        A dictionary    
    '''
    return dict(zip(rvector.names,map(list,list(rvector))))

def get_qc_class_and_string(pd_df):
    # Following ATBD from KNMI ECA&D
    n_rejects = (pd_df['p.value'] <= 0.01).sum()
    qc_classes = {
        1 : 'Useful',
        2 : 'Doubtful',
        3 : 'Suspect'
    }
    if n_rejects<=1:
        qc_class = 1
    elif n_rejects==2:
        qc_class = 2
    elif n_rejects>2:
        qc_class = 3
    return qc_class,qc_classes[qc_class]


class TrendLims1D:
    def __init__(self,name,verbose=True,params={'alpha' : 0.05, 'trend_magnitude' : None, 'n_mc' : 20000}):
        self.name = name
        self.verbose = verbose
        self.params = params
        # Initialize some empty data
        self.stats_summary = {}
        self.logbook = []
        self.fitted = {}
    def create_artificial(self,*kwargs):
        self.add_to_logbook("Created an artificial timeseries")
        # Default parameters for creating random test data
        pars_artificial = dict(periods = 100,freq = 'Y',mean = 1,trend_magnitude = 0.03,noise_magnitude = 0.05,jump_magnitude = 0.05,jump_start = 20,jump_length = 20)
        # Update the default parameters with custom parameters
        pars_artificial.update(*kwargs)

        # Create 100 years of artificial data
        timeline = pd.date_range('1/1/1901', periods=pars_artificial['periods'], freq='Y')
        # Convert to time in floats, counting in decades from starting time
        timediff = timeline-pd.datetime(1901,1,1)
        decadal_time = timediff.total_seconds()/(3600.*24.*365.25*10) # Convert to units decades

        jump_axis = np.zeros_like(decadal_time)
        jump_axis[pars_artificial['jump_start']:pars_artificial['jump_start']+pars_artificial['jump_length']] = 1
        noise_axis = np.random.uniform(-1,1,len(decadal_time))
        
        artificial_data = pars_artificial['mean'] + pars_artificial['trend_magnitude']*decadal_time + pars_artificial['jump_magnitude']*jump_axis + pars_artificial['noise_magnitude']*noise_axis

        self.data_ts = pd.Series(data=artificial_data,index=timeline)
        self.data_ts_copy = copy.deepcopy(self.data_ts)

    def load_file(self,filename):
        self.add_to_logbook("Loaded datafile {0}".format(filename))
        # Load the data
        self.data_cube = iris.load_cube(filename)
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

    def breakpoint_recipe_values(self,resamplefreq='Y'):
        '''
        # This function applies 4 homogeneity tests to yearly (by default) means of the data, similar to most climate variables in the ATBD from ECA&D
        '''
        breakpoint_input = self.data_ts.resample(resamplefreq).mean()
        self.breakpoints = self.test_breakpoint(breakpoint_input,['snh','pet','bhr','von'])
        self.qc_status_int,self.qc_status=get_qc_class_and_string(self.breakpoints)
        print("This timeseries is marked as {0} (value based)".format(self.qc_status))

    def breakpoint_recipe_differences(self,resamplefreq='Y'):
        '''
        # This function applies 4 homogeneity tests to yearly (by default) means of day-to-day absolute differences of the data, similar to one of the climate variables in the ATBD from ECA&D
        '''
        absdiff = self.data_ts.diff().abs()
        absdiff = absdiff[slice(1,None)] # Cut the first value since it is NaN
        self.differences = absdiff.resample(resamplefreq).mean()
        self.breakpoints = self.test_breakpoint(self.differences,['snh','pet','bhr','von'])
        self.qc_status_int,self.qc_status=get_qc_class_and_string(self.breakpoints)
        print("This timeseries is marked as {0} (variance based)".format(self.qc_status))

    def test_breakpoint(self,input_ts,testnames):
        values = input_ts.values
        values_as_rvector = robjects.FloatVector(values)
        rtrendpackage = rpackages.importr('trend') # This should be moved somewhere else, when doing many calls to this function
        breakpoint_results = []
        if 'snh' in testnames:
            test_result = rtrendpackage.snh_test(values_as_rvector,m=self.params['n_mc'])
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'pet' in testnames:
            test_result = rtrendpackage.pettitt_test(values_as_rvector)
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'bhr' in testnames:
            test_result = rtrendpackage.br_test(values_as_rvector,m=self.params['n_mc'])
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'bu' in testnames:
            test_result = rtrendpackage.bu_test(values_as_rvector,m=self.params['n_mc'])
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'von' in testnames:
            test_result = rtrendpackage.bartels_test(values_as_rvector)
            von_result = rvector_to_pydict(test_result)
            von_result['estimate'] = [None] # The list mimics the structure of the R output from the tests.
            von_result['estimate_formatted'] = [None] # The list mimics the structure of the R output from the tests.
            breakpoint_results.append(von_result)
        # Restructure test results to pandas dataframe
        colnames = ['method','estimate','estimate_formatted','p.value']
        rows_list = []
        for breakpoint_row in breakpoint_results:
            dict1 = OrderedDict()
            dict1.update({colname : breakpoint_row[colname][0] for colname in colnames})
            rows_list.append(dict1)
        # Create a dataframe to save breakpoint test results
        df_breakpoints = pd.DataFrame(rows_list)
        # Add additional column, with formatted breakpoints
        #df_breakpoints['estimate_formatted'] = input_ts.index[df_breakpoints['estimate']]
        return df_breakpoints

    def __get_breakpoints__(self,data_ts,resamplefreq=None):
        # By default on Monthly data
        data_ts = data_ts.resample(resamplefreq).mean()
        datavalues = data_ts.values
        data_ts_as_rvector = robjects.FloatVector(datavalues)
        rtrendpackage = rpackages.importr('trend')
        snh_result = rtrendpackage.snh_test(data_ts_as_rvector,m=self.params['n_mc'])
        snh_result = rvector_to_pydict(snh_result)
        pettitt_result = rtrendpackage.pettitt_test(data_ts_as_rvector)
        pettitt_result = rvector_to_pydict(pettitt_result)
        br_result = rtrendpackage.br_test(data_ts_as_rvector,m=self.params['n_mc'])
        br_result = rvector_to_pydict(br_result)
        bu_result = rtrendpackage.bu_test(data_ts_as_rvector,m=self.params['n_mc'])
        bu_result = rvector_to_pydict(bu_result)
        breakpoint_results = [snh_result,pettitt_result,br_result,bu_result]

        colnames = ['method','estimate','p.value']
        rows_list = []
        for test_result in breakpoint_results:
            dict1 = OrderedDict()
            dict1.update({colname : test_result[colname][0] for colname in colnames})
            rows_list.append(dict1)
        # Create a dataframe to save breakpoint test results
        df_breakpoints = pd.DataFrame(rows_list)
        # Add additional column, with formatted breakpoints
        df_breakpoints['estimate_formatted'] = data_ts.index[df_breakpoints['estimate']]
        return df_breakpoints

    def evaluate_breakpoints(self):
        # Check for each PD dataframe n_rejects
        qc_class_values,qc_string_values = get_qc_class_and_string(self.breakpoints_values)
        qc_class_absdiffvalues,qc_string_absdiffvalues = get_qc_class_and_string(self.breakpoints_absdiffvalues)

        # Print for each PD dataframe n_rejects
        print("Values: {0}".format(qc_string_values))
        print("Absolute differences: {0}".format(qc_string_absdiffvalues))

    def plot(self,label=None,mode='values'):
        self.add_to_logbook("Creating a plot with label: {0}".format(label))
        if mode=='values':
            self.data_ts.plot(label=label)
        elif mode=='differences':
            self.differences.plot(label=label)
        plt.legend()

    def plot_breakpoints(self):
        for bp_df,bp_color,shortname in zip([self.breakpoints_values,self.breakpoints_absdiffvalues],['r','b'],['on values','on absdiff']):
            xvals = bp_df['estimate_formatted'].values
            pvals = bp_df['p.value'].values
            testnames = [breakpoint_test_names[testname] for testname in bp_df['method']]
            y_base = plt.ylim()[0]
            y_range = np.abs(plt.ylim()[1]-plt.ylim()[0])
            y_offset = 0.02*y_range # 2% of the y axis range
            
            yvals = [y_base+i*y_offset for i in range(len(xvals))]
            markers = ['x','+','.','*']
            for xval,yval,marker,testname,pval in zip(xvals,yvals,markers,testnames,pvals):
                testname = testname+' '+shortname
                if pval < 0.01:
                    plt.scatter(xval,yval,marker=marker,label=testname,color=bp_color)
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),
                      ncol=3, fancybox=True, shadow=True)

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

    def trend_theilsen(self):
        from scipy.stats.mstats import theilslopes

        xaxis = self.data_ts.index.to_julian_date().values
        yaxis = self.data_ts.values

        theilsen_result = theilslopes(yaxis,x=xaxis,alpha=self.params['alpha'])
        slope,intercept,slope_low,slope_up = theilsen_result
        self.fitted['trend_theilsen'] = xaxis*slope+intercept
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

    def remove_trend(self,fit_name=None):
        self.data_ts = self.data_ts-self.fitted[fit_name]

    def do_residual_analysis(self,fit_name=None):
        '''
        # The layout of this plot follows an example by Christoph Frei in his course 'Analysis of Weather
        # and Climate Data' at ETH Zurich
        # 
        # If no fit_name is provided, the residues are calculated as the timeseries minus the mean value of the timeseries.
        '''
        if fit_name==None:
            res = self.data_ts.values-self.data_ts.values.mean()
            tukey_anscombe_xaxis = self.data_ts.index
        else:
            res = self.data_ts.values-self.fitted[fit_name]
            tukey_anscombe_xaxis = self.fitted[fit_name]

        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,squeeze=False,constrained_layout=True,figsize=(8,6))

        # ax1 QQ plot
        sm.qqplot(res/np.std(res),line='45',ax=ax1)
        ax1.set_title('Normal Q-Q Plot')
        # ax2 Histogram
        ax2.hist(res)
        ax2.set_title('Residual histogram')
        ax2.set_ylabel('# of samples')
        # testing homoscedasticity
        ax3.scatter(tukey_anscombe_xaxis,res)
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
        
