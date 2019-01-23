from c3s_511_trends import TrendLims1D
import iris

# Load test input data
filename = '/net/exo/landclim/PROJECTS/C3S/datadir/obsdir/OBS_ERA5-SM-LAYER1_ground_L3_T2Ds_sm_200001-201712.nc'
var_name = 'sm'
variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var_name))
dat3d = iris.load(filename, constraints=variable_constraint)[0]
# Downsampling
dat3d = dat3d[::10,::10,::10]
print(dat3d.shape)


import itertools as it
import datetime
from pathos.multiprocessing import ProcessingPool
import numpy as np


def custom_recipe(mydat):
    mydat.resample('Y')
    mydat.do_trends('mk')
    #mydat.remove_trend(fit_name='polynomial')
    mydat.breakpoint_recipe_values()
    # Collect all data
    data_channel = {
    'mk_pvalue' : mydat.trends['pvalue'].values[0],
    'qc_status_int' : mydat.qc_status_int,
    'algorithm_succes' : int(1),
    }
    return data_channel

def wrapperfunction(wrapdict):
        print(str(datetime.datetime.now())+": starting new task")
        #data = wrapdict['data']
        indices = wrapdict['indices']
        i,j = indices
        dat = TrendLims1D('test')
        dat.initialize_through_realization_of_cube(dat3d,i,j)
        # Data channel default
        data_channel_default = {
                'indices' : indices,
                'mk_pvalue' : np.nan,
                'qc_status_int' : np.nan,
                'algorithm_succes' : int(0),
        }
        if True:
            data_channel = custom_recipe(dat)
            data_channel_default.update(data_channel)
        #except BaseException as e:
        #    print("Gridpoint (i,j)={0},{1} returned error: {2}".format(i,j,str(e)))
        return data_channel_default


# Determine dimensions of 2D grid
shape2d = dat3d.shape[1:]
# Create tuples containing all the combinations of indices
indexlist = list(it.product(range(shape2d[0]),range(shape2d[1])))
# Dictionary to save the indices
wrapdictlist = []
for indices in indexlist:
    d = {}
    d['indices'] = indices
    wrapdictlist.append(d)

pool = ProcessingPool()
pool.ncpus = 6
output_list = pool.map(wrapperfunction,wrapdictlist)

# Now create arrays for data storage
qc_status_int = np.zeros(shape2d,dtype='float')
mk_pvalue = np.zeros(shape2d,dtype='float')
algorithm_succes = np.zeros(shape2d,dtype='int')

# Now start unpacking the results     
for d in output_list:
    i,j = d['indices']
    qc_status_int[i,j] = d['qc_status_int']
    mk_pvalue[i,j] = d['mk_pvalue']
    algorithm_succes[i,j] = d['algorithm_succes']


plt.figure(figsize=(20,15))
plt.imshow(algorithm_succes)
