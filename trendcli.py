#!/usr/bin/env python
"""Command line interface to the trend framework"""
import xarray as xr
import dask
import logging
import esmvalcore.preprocessor as pp
import argparse
import os
from trends_core3d import *
import sys

import earthreader
import matplotlib.pyplot as plt
from trendplotting import trend_mapplot

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))
logger.setLevel('DEBUG')

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--pattern',
                    '-p',
                    required=True,
                    help='pattern of input files')
parser.add_argument('--varname',
                    '-v',
                    required=True,
                    help='variable to select')
parser.add_argument('--startdate',
                    '-s',
                    default=None,
                    help='starting date as YYYY-MM')
parser.add_argument('--enddate',
                    '-e',
                    default=None,
                    help='end date as YYYY-MM')
parser.add_argument('--downsample',
                    '-d',
                    action='store_true',
                    help='Downsample horizontal resolution (one out of 10 gridpoints; for testing)')
parser.add_argument('--period',
                   '-t',
                    default='year',
                    help='Time period (year [default]; DJF; MAM; JJA; SON)')
parser.add_argument('--datasetname',
                    '-n',
                    default=None,
                    help='Name of the output file.')
parser.add_argument('--outputdir',
                    '-o',
                    required=True,
                    help='Name of the output directory (will be created).')
parser.add_argument('--statistics',
                    '-a',
                    default='mean',
                    help='Statistics to be applied over defined period.')
parser.add_argument('--alpha',
                   default=0.05,
                   type=float,
                   help='significance level for mann-kendall test')


args = parser.parse_args()

# Make output directory
datadir = f'{args.outputdir}/'
os.makedirs(datadir, exist_ok=True)


logger.info("Opening dataset.")

dask.config.set(**{'array.slicing.split_large_chunks': False})
data = earthreader.read(args.pattern, args.varname, method='xarray')
data_units = data.units

logger.info("Selecting time period.")
data = data.loc[dict(time=slice(args.startdate, args.enddate))]
if args.downsample:
    logger.info(f"Downsampling spatial resolution (1 out of 10 gridpoints).")
    data = data[:, ::10, ::10]

if args.period in ['DJF', 'MAM', 'JJA', 'SON']:
    logger.info(f"Statistics over season {args.period}")
    data = data.where(data['time.season'] == args.period)
    data = data.rolling(min_periods=3, center=True, time=3).mean()
    data = data.groupby('time.year').mean('time')
    # Rename year to time, needed for downstream functionality
    data = data.rename({'year' : 'time'})
    data.compute()
else:
    logger.info(f"Statistics over full year")
    data = data.groupby('time.year').mean('time')
    # Rename year to time, needed for downstream functionality
    data = data.rename({'year' : 'time'})
    data.compute()

logger.info(f"Calculating statistics '{args.statistics}' over period")


logger.info('Calculating theil-sen trend')
theilsen = theilsen_trend(data)
theilsen.attrs['units'] = data_units + ' year-1'
logger.info('Calculating Mann Kendall test')
mk = mannkendall(data, alpha=args.alpha)

# Mask theilsen with mk test results
theilsen_masked = theilsen.where(mk!=0, np.nan)
theilsen_masked.name += '_masked'
result = xr.merge([theilsen_masked, theilsen, mk])


corename = f"{args.varname}_theilsenmk_{args.period}{args.statistics}_{args.startdate}-{args.enddate}_alpha{str(args.alpha).replace('.','-')}"
namelist = [corename]
if args.datasetname:
    namelist.insert(0, args.datasetname)
if args.downsample:
    namelist.append('downsampled')
outname = '_'.join(namelist)+'.nc'

datasavename = os.path.join(datadir, outname)
logger.info(f'Saving: {datasavename}')
result.to_netcdf(datasavename)