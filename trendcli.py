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
                    action='store_true',
                    help='Downsample horizontal resolution (one out of 10 gridpoints; for testing)')
parser.add_argument('--period',
                   '-t',
                    default='year',
                    help='Time period (year; DJF; MAM; JJA; SON)')
parser.add_argument('--datasetname',
                    '-n',
                    default=None,
                    help='Name of the output file.')
parser.add_argument('--statistics',
                    '-a',
                    default='mean',
                    help='Statistics to be applied over defined period.')
parser.add_argument('--alpha',
                   default=0.05,
                   type=float,
                   help='significance level for mann-kendall test')

args = parser.parse_args()



logger.info("Opening dataset.")
data = xr.open_mfdataset(args.pattern, combine='by_coords')
data = data[args.varname]
logger.info("Selecting time period.")
data = data.loc[dict(time=slice(args.startdate, args.enddate))]
if args.downsample:
    logger.info(f"Downsampling spatial resolution (1 out of 10 gridpoints).")
    data = data[:, ::10, ::10]

# Now convert to Iris for using esmvalcore preprocessor functions
logger.debug("Converting to Iris cube")
cube = data.to_iris()

if args.period in ['DJF', 'MAM', 'JJA', 'SON']:
    logger.info(f"Statistics over season {args.period}")
    cube = pp.extract_season(cube, args.period)
else:
    logger.info(f"Statistics over full year")

logger.info(f"Calculating statistics '{args.statistics}' over period")
cube = pp.annual_statistics(cube, operator=args.statistics)

# Back to xarray
data = xr.DataArray.from_iris(cube)
logger.info('Calculating theil-sen trend')
theilsen = theilsen_trend(data)
logger.info('Calculating Mann Kendall test')
mk = mannkendall(data, alpha=args.alpha)

# Mask theilsen with mk test results
theilsen = theilsen.where(mk!=0, np.nan)
result = xr.merge([theilsen, mk])


corename = f"{args.period}{args.statistics}_{args.startdate}-{args.enddate}_alpha{args.alpha}"
namelist = [corename]
if args.datasetname:
    namelist.insert(0, args.datasetname)
if args.downsample:
    namelist.append('downsampled')
outname = '_'.join(namelist)+'.nc'

logger.info(f'Saving: {outname}')
result.to_netcdf(outname)