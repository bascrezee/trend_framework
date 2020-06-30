import numpy as np
import itertools as it
import dask
import xarray as xr
from diag1d import *

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
    # Return all nan if any nan in array
    if np.any(~np.isfinite(inputdata)):
        result = np.ones_like(inputdata)
        result[:] = np.nan
        return result

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
    if 'units' in inputdata.attrs:
        lintrend.attrs['units'] = inputdata.attrs['units'] + ' / timestep'
    else:
        lintrend.attrs['units'] = 'per timestep'
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
    if theilsentrend.name:
        theilsentrend.name += '_theilsenslope'
    else:
        theilsentrend.name = 'theilsenslope'
    if 'units' in inputdata.attrs:
        theilsentrend.attrs['units'] = inputdata.attrs['units'] + ' per timestep'
    else:
        theilsentrend.attrs['units'] = 'per timestep'
    return theilsentrend

def mannkendall(inputdata, alpha=0.05):
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
    mannkendall_da = dask.array.apply_along_axis(mannkendall1d,0,inputdata, alpha=alpha)
    template = inputdata[0:1].mean('time')
    mannkendall = template.copy()
    mannkendall.values = mannkendall_da
    if mannkendall.name:
        mannkendall.name += '_mk_test'
    else:
        mannkendall.name = 'mk_test'
    mannkendall.attrs['info'] = '-1: a negative trend; 0: no trend; 1: a positive trend'
    mannkendall.attrs['units'] = 1
    return mannkendall




# Helper functions
def validfrac(inputdata,axis=0):
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
    validfrac_da = dask.array.apply_along_axis(validfrac1d,axis,inputdata)
    template = inputdata[0:1].mean('time')
    validfrac = template.copy()
    validfrac.values = validfrac_da
    validfrac.name += '_valid_fraction'
    validfrac.attrs['units'] = 1
    return validfrac


# Calculate valid fraction of aggregated data
def validfrac_along_axis(x,axis):
    return np.apply_along_axis(lambda y: np.isfinite(y).sum()/y.size,axis,x)

# Advanced xarray helper functions
def broadcast_time(inputdata,whichtime='time.dayofyear'):
    '''
    This function broadcasts time over the full dimension of the dataset

    The function is intended to check if changes in temporal coverage could impact the trend tests
    '''
    myds = inputdata.to_dataset()
    doy,dummy = xr.broadcast(myds[whichtime],myds)
    doy_masked = doy.where(np.isfinite(inputdata.values))
    return doy_masked
