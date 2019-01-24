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
import rpy2.rinterface
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
    try:
        result = dict(zip(rvector.names,map(list,list(rvector))))
    except TypeError:
        print("Warning: setting RNULLType to FloatVector with length zero before conversion to Python dictionary")
        # probably there is a RNULLType in the object, which can not be converted
        # as a workaround, convert it to an empty FloatVector
        for i in range(len(rvector)):
            if type(rvector[i])==rpy2.rinterface.RNULLType:
                rvector[i]=robjects.FloatVector([])
        result = dict(zip(rvector.names,map(list,list(rvector))))
    return result

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
    '''  Object representing one-dimensional timeseries data and specified methods for detection of trends, breakpoints etc. 
    '''

    def __init__(self,name,verbose=True):
        '''Initialization of the object with the option to set some parameters that are used globally

        Parameters
        ----------
        name : str
            Name of the object
        '''
        self.name = name
        self.verbose = verbose
        # Initialize some empty data
        self.stats_summary = {}
        self.logbook = []
        self.fitted = {}

    def create_artificial(self,*kwargs):
        '''Function to create artificial data

        Example usage: create_artificial(periods = 100,freq = 'Y',mean = 1,trend_magnitude = 0.03,noise_magnitude = 0.05,jump_magnitude = 0.05,jump_start = 20,jump_length = 20)
        '''
        self.__add_to_logbook__("Created an artificial timeseries")
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

    def load_file(self,filename,var_name=None):
        ''' Loading a netCDF file. If the file has more than one variable, the variable name needs to be given 
        
        Parameters
        ----------
            filename : str
                the filename
            var_name : str, optional
                the variable name (only needed when more than one variable is present on the file)
        '''
        self.__add_to_logbook__("Loaded datafile {0}".format(filename))
        # Load the data
        if var_name is not None:
            variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var_name))
            self.data_cube = iris.load(filename, constraints=variable_constraint)[0]
        else:
            self.data_cube = iris.load_cube(filename)
        self.data_cube = iris.util.squeeze(self.data_cube)
        self.data_ts = iris.pandas.as_series(self.data_cube)
        self.data_ts_copy = copy.deepcopy(self.data_ts)
        self.data_cube_copy = copy.deepcopy(self.data_cube)

        # Do a check on missing values
        self.__check_missingvalues__()

    def initialize_through_realization_of_cube(self,cube3d,i,j):
        ''' This function takes a certain gridpoint from a lazy cube and realizes this gridpoint '''
        self.data_cube = iris.util.squeeze(cube3d[:,i,j])
        self.data_ts = iris.pandas.as_series(self.data_cube)
        # Do not add a copy in this case

    def __check_missingvalues__(self):
        if self.data_ts.isnull().any():
            print("Warning: the data contains missing values")

    def reset(self):
        ''' Restore the timeseries data to the state right after reading in the file or having created the artificial data
        '''
        self.data_cube = self.data_cube_copy
        self.data_ts = self.data_ts_copy
        self.logbook = self.logbook[1:] # Only preserve first line, with info on loading

    def subset(self,timeslice):
        ''' Subset over time.

        Parameters
        ----------
            timeslice : slice
                a slice object representing the wanted timeslice
        
        Example usage: 
            mydat.subset(slice('1950-01-01','2017-12-31'))

        '''
        self.data_ts = self.data_ts[timeslice]
        self.__add_to_logbook__('Subsetted to timeperiod {0}-{1}'.format(timeslice.start,timeslice.stop))
        # Do a check on missing values
        self.__check_missingvalues__()

    def extract_season(self,season):
        ''' Extract a certain season. 

        NB: Apply this function right after loading the data

        Parameters
        ----------
            season : str
                either 'djf','mam','jja','son'
        '''
        import iris.coord_categorisation
        clim_seasons = ['djf','mam','jja','son']
        if season not in clim_seasons:
            print("Season should be one of the following: ",clim_seasons)
            raise ValueError
        try:
            iris.coord_categorisation.add_season(self.data_cube, 'time', name='clim_season')
        except ValueError: # Probably the coordinate has been added already, therefore pass. 
            pass
        constraint = iris.Constraint(clim_season=season)
        self.data_cube = self.data_cube.extract(constraint)
        # Update timeseries as well
        self.data_ts = iris.pandas.as_series(self.data_cube)
        self.__add_to_logbook__('Extracted season {0}'.format(season))

    def extract_month(self,month):
        ''' Extract a certain month.

        NB: Apply this function right after loading the data

        Parameters
        ----------
            month : str
                either 'Jan','Feb','Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'
        '''
        import iris.coord_categorisation
        clim_months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
        if month not in clim_months:
            print("Month should be one of the following: ",clim_months)
            raise ValueError
        try:
            iris.coord_categorisation.add_month(self.data_cube, 'time', name='months')
        except ValueError: # Probably the coordinate has been added already, therefore pass. 
            pass
        constraint = iris.Constraint(months=month)
        self.data_cube = self.data_cube.extract(constraint)
        # Update timeseries as well
        self.data_ts = iris.pandas.as_series(self.data_cube)
        self.__add_to_logbook__('Extracted month {0}'.format(month))

    def decompose_seasonal(self,inplace=False,freq=None):
        ''' Decompose the data into a seasonal cycle, a trend, and a residual. If used with inplace=True, it deseasonalizes the data.

        Parameters
        ----------
            inplace : bool
                If False, the result is returned. If True, the object is deseasonalized
        '''
        results = seasonal_decompose(self.data_ts,model='additive',freq=freq)
        signal = results.trend
        seasonal = results.seasonal
        residue = results.resid
        # Note that the seasonal fit, is a sum of the seasonal sign and the dataset mean
        self.fitted['seasonal'] = results.seasonal 
        if inplace:
            self.data_ts = self.data_ts-(self.fitted['seasonal'])
            self.__add_to_logbook__('Deseasonalized the data')
        else:
            return results

    def handle_missingdata(self,strategy='dropna'):
        ''' This function takes care of the missing data. 
        
        Parameters
        ----------
        strategy : str
            The strategy for handling missing data (either dropna or fillgaps)
        '''
        if strategy=='dropna':
            self.data_ts = self.data_ts.dropna()
        elif strategy=='fillgaps':
            raise NotImplementedError

    def resample(self,resamplefreq):
        ''' Resample the data to the target frequency

        Parameters
        ----------
        resamplefreq : str
            The frequency to which to resample the data

        Example usage: 
            mydat.resample('Y') # Y=Yearly, M=Monthly, D=Daily
        '''
        self.data_ts = self.data_ts.resample(resamplefreq).mean()
        # Do a check on missing values
        self.__check_missingvalues__()
        self.__add_to_logbook__('Resampled to {0} frequency'.format(resamplefreq))

    def breakpoint_recipe_values(self,resamplefreq='Y',rpackagename='trend'):
        ''' Apply four homogeneity tests to yearly means of the data, similar to most climate variables in the ATBD from ECA&D

        The tests either accept or reject the null hypothesis (no breakpoints present):
        
        If n_rejects <= 1        --> data marked as usefull
        Else if n_rejects == 2   --> data marked as doubtfull
        Else if n_rejects >2     --> data marked as suspect

        Parameters
        ----------
        resamplefreq : str
            The frequency to which to resample the data before applying the breakpoint detection
        rpackagename : str
            The name of the R-package to use for the breakpoint detection, either 'trend' (default) or 'iki.dataclim'

        References
        ----------
        rpackagename='trend':
            https://cran.r-project.org/web/packages/trend/index.html
        rpackagename='iki.dataclim':
            https://cran.r-project.org/web/packages/iki.dataclim/index.html
        '''
        breakpoint_input = self.data_ts.resample(resamplefreq).mean()
        if rpackagename=='trend':
            self.breakpoints = self.__test_breakpoint_rtrend__(breakpoint_input,['snh','pet','bhr','von'])
            self.qc_status_int,self.qc_status=get_qc_class_and_string(self.breakpoints)
            print("This timeseries is marked as {0} (value based)".format(self.qc_status))
        elif rpackagename=='iki.dataclim':
            self.breakpoints = self.__test_breakpoint_rikidataclim__(breakpoint_input,['snh','pet','bhr','von'])

    def breakpoint_recipe_differences(self,resamplefreq='Y',rpackagename='trend'):
        ''' Apply four homogeneity tests to yearly means of day-to-day absolute differences of the data, similar to most climate variables in the ATBD from ECA&D

        The tests either accept or reject the null hypothesis (no breakpoints present):
        
        If n_rejects <= 1        --> data marked as usefull
        Else if n_rejects == 2   --> data marked as doubtfull
        Else if n_rejects >2     --> data marked as suspect

        Parameters
        ----------
        resamplefreq : str
            The frequency to which to resample the data before applying the breakpoint detection
        rpackagename : str
            The name of the R-package to use for the breakpoint detection, either 'trend' (default) or 'iki.dataclim'

        References
        ----------
        rpackagename='trend':
            https://cran.r-project.org/web/packages/trend/index.html
        rpackagename='iki.dataclim':
            https://cran.r-project.org/web/packages/iki.dataclim/index.html
        '''
        absdiff = self.data_ts.diff().abs()
        absdiff = absdiff[slice(1,None)] # Cut the first value since it is NaN
        self.differences = absdiff.resample(resamplefreq).mean()
        if rpackagename=='trend':
            self.breakpoints = self.__test_breakpoint_rtrend__(self.differences,['snh','pet','bhr','von'])
            self.qc_status_int,self.qc_status=get_qc_class_and_string(self.breakpoints)
            print("This timeseries is marked as {0} (variance based)".format(self.qc_status))
        elif rpackagename=='iki.dataclim':
            self.breakpoints = self.__test_breakpoint_rikidataclim__(self.differences,['snh','pet','bhr','von'])

    def __test_breakpoint_rtrend__(self,input_ts,testnames,n_mc=20000):
        values = input_ts.values
        values_as_rvector = robjects.FloatVector(values)
        rtrendpackage = rpackages.importr('trend') # This should be moved somewhere else, when doing many calls to this function
        breakpoint_results = []
        if 'snh' in testnames:
            test_result = rtrendpackage.snh_test(values_as_rvector,m=n_mc)
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'pet' in testnames:
            test_result = rtrendpackage.pettitt_test(values_as_rvector)
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'bhr' in testnames:
            test_result = rtrendpackage.br_test(values_as_rvector,m=n_mc)
            test_result = rvector_to_pydict(test_result)
            test_result['estimate_formatted'] = [input_ts.index[test_result['estimate'][0]]]
            breakpoint_results.append(test_result)
        if 'bu' in testnames:
            test_result = rtrendpackage.bu_test(values_as_rvector,m=n_mc)
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

    def __test_breakpoint_rikidataclim__(self,input_ts,testnames):
        values = input_ts.values
        values_as_rvector = robjects.FloatVector(values)
        rtrendpackage = rpackages.importr('iki.dataclim') # This should be moved somewhere else, when doing many calls to this function
        breakpoint_results = []
        if 'snh' in testnames:
            test_result = rtrendpackage.SNHtest(values_as_rvector)
            test_result = rvector_to_pydict(test_result)
            test_result['method']='SNH Test'
            breakpoint_results.append(test_result)
        if 'pet' in testnames:
            test_result = rtrendpackage.PETtest(values_as_rvector)
            test_result = rvector_to_pydict(test_result)
            test_result['method']='Petitt Test'
            breakpoint_results.append(test_result)
        if 'bhr' in testnames:
            test_result = rtrendpackage.BHRtest(values_as_rvector)
            test_result = rvector_to_pydict(test_result)
            test_result['method']='Buishand Range Test'
            breakpoint_results.append(test_result)
        if 'von' in testnames:
            test_result = rtrendpackage.VONtest(values_as_rvector)
            von_result = rvector_to_pydict(test_result)
            von_result['method']='Von Neumann Ratio test'
#           von_result['estimate'] = [None] # The list mimics the structure of the R output from the tests.
#            von_result['estimate_formatted'] = [None] # The list mimics the structure of the R output from the tests.
            breakpoint_results.append(von_result)
        # Restructure test results to pandas dataframe
        colnames = ['method','statistic','breakpoint','significance']
        rows_list = []
        for breakpoint_row in breakpoint_results:
            dict1 = OrderedDict()
            dict1.update({colname : breakpoint_row[colname] for colname in colnames})
            rows_list.append(dict1)
        # Create a dataframe to save breakpoint test results
        df_breakpoints = pd.DataFrame(rows_list)
        # Add additional column, with formatted breakpoints
        #df_breakpoints['estimate_formatted'] = input_ts.index[df_breakpoints['estimate']]
        return df_breakpoints

    def __evaluate_breakpoints__(self):
        # Check for each PD dataframe n_rejects
        qc_class_values,qc_string_values = get_qc_class_and_string(self.breakpoints_values)
        qc_class_absdiffvalues,qc_string_absdiffvalues = get_qc_class_and_string(self.breakpoints_absdiffvalues)

        # Print for each PD dataframe n_rejects
        print("Values: {0}".format(qc_string_values))
        print("Absolute differences: {0}".format(qc_string_absdiffvalues))

    def plot(self,label=None,mode='values',marker='.',ms=4):
        self.__add_to_logbook__("Creating a plot with label: {0}".format(label))
        if mode=='values':
            self.data_ts.plot(label=label,marker=marker,ms=ms)
        elif mode=='differences':
            self.differences.plot(label=label,marker=marker,ms=ms)
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

    def trend_mktest(self,alpha=0.05):
        mk_input = self.data_ts.values
        mk_input = mk_input[~np.isnan(mk_input)]
        trend,h,p,z = self.__mk_test__(mk_input,alpha=alpha)
        # Store results in stats_summary
        results = {}
        results['sign']=trend
        results['slope']=None
        results['h']=h
        results['pvalue']=p
        results['z']=z
        self.__add_to_logbook__('Calculated Mann-Kendall test')
        results['method']='mk'
        return results

    def trend_linear(self,alpha=0.05):
        xaxis = self.data_ts.index.to_julian_date()
        yaxis = self.data_ts.values
        results = linregress(xaxis,yaxis)
        self.fitted['linear'] = xaxis*results.slope+results.intercept
        results_dict = {}
        results_dict['slope']=results.slope*(365.25*10) # From /day to /decade
        if results.pvalue <= alpha:
            results_dict['sign'] = int(np.sign(results.slope))
        else:
            results_dict['sign'] = int(0)
        results_dict['pvalue']=results.pvalue
        results_dict['stderr']=results.stderr*(365.25*10)
        results_dict['slope_low'] = results_dict['slope']-results_dict['stderr']
        results_dict['slope_up'] = results_dict['slope']+results_dict['stderr']
        self.__add_to_logbook__('Calculated linear trend test')
        results_dict['method']='linear'
        return results_dict

    def trend_theilsen(self,alpha=0.05):
        from scipy.stats.mstats import theilslopes

        xaxis = self.data_ts.index.to_julian_date().values
        yaxis = self.data_ts.values

        theilsen_result = theilslopes(yaxis,x=xaxis,alpha=alpha)
        slope,intercept,slope_low,slope_up = theilsen_result
        self.fitted['theilsen'] = xaxis*slope+intercept
        assert(slope_low <= slope <= slope_up) # Just to be safe, check this
        slope_sign = np.sign(slope)
        if not slope_low < 0.0 < slope_up:
            sign = int(np.sign(slope))
        else:
            sign = int(0)
        
        results_dict = {}
        results_dict['sign'] = sign
        results_dict['slope'] = slope*(365.25*10) # From /day to /decade
        results_dict['slope_low'] = slope_low*(365.25*10) # From /day to /decade
        results_dict['slope_up'] = slope_up*(365.25*10) # From /day to /decade
        # Add a sign to Theilsen


        self.__add_to_logbook__('Calculated Theil-Sen slope')
        results_dict['method']='theilsen'
        results_dict['pvalue']=None # Trend Theilsen has no pvalue
        return results_dict
        
    def do_trends(self,trend_names=['mk','linear','theilsen'],alpha=0.05):
        ''' This function calls the different trend tests.

        Parameters
        ----------
        trend_names : str
            The names of the trend tests that need to be calculated

        The units of the trends are 'variable unit per decade' and are saved in trends

        Returns
        -------
            A dictionary containing:
                trend_magnitude : same as above
        '''
        trend_results = []
        if 'mk' in trend_names:
            mk_result = self.trend_mktest(alpha=alpha)
            trend_results.append(mk_result)
        if 'linear' in trend_names:
            linear_result = self.trend_linear(alpha=alpha)
            trend_results.append(linear_result)
        if 'theilsen' in trend_names:
            theilsen_result = self.trend_theilsen(alpha=alpha)
            trend_results.append(theilsen_result)
        # Restructure test results to pandas dataframe
        colnames = ['method','sign','slope','pvalue']
        rows_list = []
        for trend_row in trend_results:
            dict1 = OrderedDict()
            dict1.update({colname : trend_row[colname] for colname in colnames})
            rows_list.append(dict1)
        # Create a dataframe to save trend test results
        df_trends = pd.DataFrame(rows_list)
        # Rename column(s)
        # Quick fix. Look into issue
        try:
            # Calculate additionally the slope in percentage
            df_trends['slope [percent/decade]'] = df_trends['slope']/self.data_cube.collapsed('time',iris.analysis.MEAN).data
            # First reorder the columns
            columns_right_order = ['method','sign','slope','slope [percent/decade]','pvalue']
            df_trends = df_trends[columns_right_order]
            # Then rename a certain column
            slope_name = 'slope [{0}/decade]'.format(str(self.data_cube.units))
            df_trends.rename(columns={'slope' : slope_name},inplace=True)
        except AttributeError:
            pass
        self.trends = df_trends

    def remove_trend(self,fit_name=None):
        '''
        Remove a trend from the data. In case a linear trend on the data is expected, one 
        can use 'theilsen' or 'linear'. 
        For non-linear trends one can choose between either 'polynomial' or 'differences' where 
        the latter should be chosen in case of a fully stochastic trend.
        '''
        if fit_name in ['linear','theilsen']:
            self.data_ts = self.data_ts-self.fitted[fit_name]
        elif fit_name=='polynomial':
            xaxis = self.data_ts.index.to_julian_date().values
            yaxis = self.data_ts.values
            d,c,b,a = np.polyfit(xaxis,yaxis,3)
            y_fit = a+b*xaxis+c*xaxis**2+d*xaxis**3
            self.data_ts = self.data_ts-y_fit
        else: 
            print("Invalid fit_name: ",fit_name)
            raise ValueError

    def weatherhead_framework(self,trend_name='theilsen',trend_magnitude=None):
        ''' This framework follows  Weatherhead et al. 1998 [1] for estimating the amount of years needed to detect a trend of certain magnitude with a probability of 90%.

        Data is initially resampled to yearly means.

        Parameters
        ----------
        trend_name : str
            The name of the trend method used for detrending the data
        trend_magnitude : float
            The magnitude of the trend in [data units/decade] in which one would be interested

        Returns
        -------
            A dictionary containing:
                trend_magnitude : same as above
                std_res : the standard deviation of the residuals
                acf_res : the 1-lag autocorrelation
                n_star : the estimated number of years
        
        References
        -----------
        [1] Weatherhead, Betsy & C. Reinsel, Gregory & C. Tiao, George & Meng, Xiao-Li & Choi, Dongseok & Cheang, Wai-Kwong & Keller, Teddie & DeLuisi, John & Wuebbles, Donald & Kerr, J & J. Miller, Alvin & Oltmans, Samuel. (1998). Factors affecting the detection of trends: Statistical considerations and applications to environmental data. Journal of Geophysical Research. 1031. 17149-17162. 10.1029/98JD00995. 
        '''
        self.data_ts = self.data_ts.resample('Y').mean()
        self.do_trends(trend_names=[trend_name])
        self.remove_trend(fit_name=trend_name)
        # Calculate according to Weatherhead et al. 
        std_res = np.std(self.data_ts.values)
        acf_res = acf(self.data_ts.values,nlags=1)[-1]
        trend_magnitude_per_year = trend_magnitude/10.
        n_star = ((3.3*std_res/np.abs(trend_magnitude_per_year))*((1+acf_res)/(1-acf_res))**.5)**(2./3)
        weatherhead_dict = OrderedDict({'trend_magnitude [/decade]' : trend_magnitude, 'std_res' : std_res,'acf_res' : acf_res, 'n_star' : n_star},index=[0])

        # Restructure test results to pandas dataframe
        # Create a dataframe to save weatherhead framework results
        self.weatherhead = pd.DataFrame(weatherhead_dict)
        return self.weatherhead

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

    def __mk_test__(self,x, alpha=0.05):
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

    def __add_to_logbook__(self,logmessage):
        if self.verbose:
            print(logmessage)
        self.logbook.append(logmessage+'\n')

    def print_logbook(self):
        print('\n'.join(self.logbook))

    def print_stats(self):
        pprint.pprint(self.stats_summary)
        










def assess_trend_consistency(trendobject):
    trend_consistency_score = 0
    all_trends = [trendobject.stats_summary[trendname]['trend'] for trendname in ['linear','mk','theilsen']]
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
    name_a,name_b='linear','theilsen'
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
        
