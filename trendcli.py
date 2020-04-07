#!/net/exo/landclim/crezees/conda/envs/jupyterlab/bin/python
"""Command line interface to the trend framework""" 
import xarray as xr
import logging
import esmvalcore.preprocessor as pp
import argparse
from trends_core3d import *
import sys

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))
logger.setLevel('DEBUG')

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--pattern',
                    '-p',
                    default=None,
                    help='pattern of input files')
parser.add_argument('--varname',
                    '-v',
                    default=None,
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
                    default=None,
                    type=int,
                    help='Downsample horizontal resolution by taking only 1 in d gridpoints')
parser.add_argument('--seltime',
                   '-t',
                    default='year',
                    help='Time aggregation (year; DJF; MAM; JJA; SON)')
parser.add_argument('--datasetname',
                    '-n',
                    default=None,
                    help='Name of the output file.')
parser.add_argument('--aggregation',
                    '-a',
                    default='mean',
                    help='Statistics to be applied.')


args = parser.parse_args()

logger.info("Opening dataset.")
data = xr.open_mfdataset(args.pattern, combine='by_coords')
data = data[args.varname]
logger.info("Selecting time period.")
data = data.loc[dict(time=slice(args.startdate, args.enddate))]
if args.downsample:
    logger.info(f"Downsampling spatial resolution (1 out of {args.downsample} gridpoints).")
    data = data[:, ::args.downsample, ::args.downsample]

# Now convert to Iris for nice preprocessor functions
logger.debug("Converting to Iris cube")
cube = data.to_iris()

if args.seltime in ['DJF', 'MAM', 'JJA', 'SON']:
    logger.info(f"Selecting season {args.seltime}")
    cube = pp.extract_season(cube, args.seltime)
else:
    logger.info(f"Statistics over full year")

logger.info(f"Calculating {args.aggregation} over period")
cube = pp.annual_statistics(cube, operator=args.aggregation)

# Back to xarray
data = xr.DataArray.from_iris(cube)
logger.info('Calculating theil-sen trend')
theilsen = theilsen_trend(data)
logger.info('Calculating Mann Kendall test')
mk = mannkendall(data)

result = xr.merge([theilsen, mk])

if args.datasetname:
    outname = f'trend_{args.datasetname}_{args.timemean}_{args.startdate}_{args.enddate}.nc'
else:
    outname = f'trend_{args.seltime}_{args.startdate}_{args.enddate}.nc'
logger.info(f'Saving: {outname}')
result.to_netcdf(outname)