import numpy as np
from scipy.stats.mstats import linregress
from scipy.stats import norm
import itertools as it
import dask
from scipy.signal import detrend
from statsmodels.tsa.stattools import acf

def weatherhead(inputdata, trend_magnitude=None):
    '''
    This function applies the Weatherhead framework (see function weatherhead1d) to a xarray DataArray
    
    Parameters
    ----------
    inputdata : xarray DataArray
        the DataArray to which the function should be applied
    trend_magnitude : float
        the trend magnitude (units: per year) in which one is interested
    
    Returns
    -------
    three xarray DataArrays:
        standard deviation ; lag-1 autocorellation ; n_star 
    '''
    std_res_da, acf_res_da, n_star_da = dask.array.apply_along_axis(weatherhead1d,0,inputdata,trend_magnitude=trend_magnitude)
    # Now put results into a DataArray
    # For std_res
    template = inputdata[0:1].mean('time')
    std_res = template.copy()
    std_res.values = std_res_da
    std_res.name += 'std of {0}'.format(template.name)
    # For acf_res
    acf_res = template.copy()
    acf_res.values = acf_res_da
    acf_res.name += 'acf of {0}'.format(template.name)
    # For n_star
    n_star = template.copy()
    n_star.values = n_star_da
    n_star.name = 'nstar of {0}'.format(template.name)
    n_star.attrs['units'] = 'years'
    return std_res,acf_res,n_star
    
def remove_trend(inputdata):
    '''
    This function detrends an xarray DataArray using a gridpoint-wise least square error linear fit (see scipy.signal.detrend). 
    
    Parameters
    ----------
    inputdata : xarray DataArray
        the DataArray that should be detrended
    
    Returns
    -------
    xarray DataArray:
        detrended dataset
    '''
    inputdata_detrended = inputdata.copy()
    inputdata_detrended.values = dask.array.apply_over_axes(detrend,inputdata,[0])
    return inputdata_detrended

def linear_trend(inputdata):
    '''
    This function calculates linear trend and p-value from an xarray DataArray
    
    Parameters
    ----------
    inputdata : xarray DataArray
    
    Returns
    -------
    xarray DataArray (two times):
        linear trend ; p-value
    '''
    lintrend_da,linpvalue_da = dask.array.apply_along_axis(lineartrend1d,0,inputdata)
    # Now put results into a DataArray
    # For the linear trend itself
    template = inputdata[0:1].mean('time')
    lintrend = template.copy()
    lintrend.values = lintrend_da
    lintrend.name += '_linslope'
    lintrend.attrs['units'] = inputdata.attrs['units'] + ' / timestep'
    # For its pvalue
    linpvalue = template.copy()
    linpvalue.values = linpvalue_da
    linpvalue.name = 'p-value of linear trend test'
    linpvalue.attrs['units'] = 1
    return lintrend, linpvalue

def theilsen_trend(inputdata):
    '''
    This function calculates theilsen slopes from an xarray DataArray
    
    Parameters
    ----------
    inputdata : xarray DataArray
    
    Returns
    -------
    xarray DataArray:
        theilsen slope
        
    see also: theilslopes1d
    '''
    theilsentrend_da = dask.array.apply_along_axis(theilslopes1d,0,inputdata)
    template = inputdata[0:1].mean('time')
    theilsentrend = template.copy()
    theilsentrend.values = theilsentrend_da
    theilsentrend.name += '_theilsenslope'
    theilsentrend.attrs['units'] = inputdata.attrs['units'] + ' / timestep'
    return theilsentrend

def mannkendall(inputdata):
    '''
    This function calculates Mann-Kendall trend test from an xarray DataArray
    
    Parameters
    ----------
    inputdata : xarray DataArray
    
    Returns
    -------
    xarray DataArray:
        Mann Kendall test outcome
        
    see also: mannkendall1d
    '''
    mannkendall_da = dask.array.apply_along_axis(mannkendall1d,0,inputdata)
    template = inputdata[0:1].mean('time')
    mannkendall = template.copy()
    mannkendall.values = mannkendall_da
    mannkendall.name += '_mk_test'
    mannkendall.attrs['info'] = '-1: a negative trend; 0: no trend; 1: a positive trend'
    mannkendall.attrs['units'] = 1
    return mannkendall

def weatherhead1d(inputdata,trend_magnitude=None):
    ''' This framework follows  Weatherhead et al. 1998 [1] for estimating the amount of years needed to detect a trend of certain magnitude with a probability of 90%.
    Data has to be provided with yearly frequency (e.g. yearly means, means for DJF, etc.)

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
    # Calculate according to Weatherhead et al.
    std_res = np.std(inputdata)
    # This is a workaround, related to the following Dask issue: https://github.com/dask/dask/pull/3742
    if len(inputdata)==1:
        acf_res = inputdata
    else:
        acf_res = acf(inputdata, nlags=1)[1]
    n_star = ((3.3 * std_res / np.abs(trend_magnitude)) * (
        (1 + acf_res) / (1 - acf_res))**.5)**(2. / 3.)
    return std_res, acf_res, n_star

def mannkendall1d(x, alpha=0.05):
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
    x = x[~np.isnan(x)] # Check this.
    
    n = len(x)

    # calculate S
    s = 0
    for k in range(n - 1):
        for j in range(k + 1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n * (n - 1) * (2 * n + 5)) / 18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(x == unique_x[i])
        var_s = (n * (n - 1) * (2 * n + 5) - np.sum(tp * (tp - 1) *
                                                    (2 * tp + 5))) / 18

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:  # s == 0:
        z = 0

    # calculate the p_value
    p = 2 * (1 - norm.cdf(abs(z)))  # two tail test
    h = abs(z) > norm.ppf(1 - alpha / 2)

    if (z < 0) and h:
        trend = -1
    elif (z > 0) and h:
        trend = 1
    else:
        trend = 0
    return trend #, h, p, z


def lineartrend1d(y,x=None, alpha=0.05):
    y = np.array(y).flatten()
    if x is None:
        x = np.arange(len(y), dtype=float)
    else:
        x = np.array(x, dtype=float).flatten()
    linoutput = linregress(x, y)
    return linoutput.slope, linoutput.pvalue

def theilslopes1d(y,x=None):
    '''
    Adapted from scipy.stats.theilslopes, leaving out calculation of confidence intervals and allowing for nan values.
    '''
    y = np.array(y).flatten()
    if x is None:
        x = np.arange(len(y), dtype=float)
    else:
        x = np.array(x, dtype=float).flatten()
    if len(x) != len(y):
        raise ValueError("Incompatible lengths ! (%s<>%s)" % (len(y), len(x)))
    # Compute sorted slopes only when deltax > 0
    deltax = x[:, np.newaxis] - x
    deltay = y[:, np.newaxis] - y
    slopes = deltay[deltax > 0] / deltax[deltax > 0]
    medslope = np.nanmedian(slopes)
    return medslope


# Helper functions
def validfrac(inputdata):
    '''
    This function calculates the valid fraction along the time axis from an xarray DataArray based on the numpy function `isfinite()`
    
    Parameters
    ----------
    inputdata : xarray DataArray
    
    Returns
    -------
    xarray DataArray:
        valid fraction
        
    see also: np.isfinite()
    '''
    validfrac_da = dask.array.apply_along_axis(validfrac1d,0,inputdata)
    template = inputdata[0:1].mean('time')
    validfrac = template.copy()
    validfrac.values = validfrac_da
    validfrac.name += '_valid_fraction'
    validfrac.attrs['units'] = 1
    return validfrac

def validfrac1d(x):
    return np.isfinite(x).sum()/x.size