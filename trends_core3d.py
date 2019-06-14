import numpy as np
from scipy.stats.mstats import linregress
from scipy.stats import norm
import itertools as it
import dask
from scipy.signal import detrend
from statsmodels.tsa.stattools import acf
        
def remove_trend(inputdata):
    inputdata_detrended = inputdata.copy()
    inputdata_detrended.values = dask.array.apply_over_axes(detrend,inputdata,[0]).compute()
    return inputdata_detrended

def linear_trend(inputdata):
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
    theilsentrend_da = dask.array.apply_along_axis(theilslopes1d,0,inputdata)
    template = inputdata[0:1].mean('time')
    theilsentrend = template.copy()
    theilsentrend.values = theilsentrend_da
    theilsentrend.name += '_theilsenslope'
    theilsentrend.attrs['units'] = inputdata.attrs['units'] + ' / timestep'
    return theilsentrend

def mannkendall(inputdata):
    mannkendall_da = dask.array.apply_along_axis(mannkendall1d,0,inputdata)
    template = inputdata[0:1].mean('time')
    mannkendall = template.copy()
    mannkendall.values = mannkendall_da
    mannkendall.name += '_mk_test'
    mannkendall.attrs['info'] = '-1: a negative trend; 0: no trend; 1: a positive trend'
    mannkendall.attrs['units'] = 1
    return mannkendall


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