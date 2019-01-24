
# coding: utf-8

# In[1]:


import sys
sys.path.append("..")
from c3s_511_trends import TrendLims1D


# # Example A: Reading in station data, resampling and plotting
# 

# In[2]:


mydat = TrendLims1D('De Bilt T2m')
mydat.load_file('../test_data/ehdb_t2m.nc')
mydat.resample('Y')
mydat.plot('T2m')


# # Example B1: Creating artificial data and breakpoint detection

# In[3]:


mydat = TrendLims1D('artificial_timeseries')
mydat.create_artificial({'periods' : 200, 'jump_start' : 46, 'jump_length': 200,'jump_magnitude' : 0.05, 'trend_magnitude' : 0.005})
mydat.do_trends()
mydat.resample('Y') # Resample to yearly resolution
mydat.plot(label='Test data')
mydat.detect_breakpoints()


# # Example B2: How breakpoint detection is affected by removal of the linear trend

# In[4]:


mydat = TrendLims1D('artificial_timeseries')
print("Creating test data with no jump")
noise_magnitude = 0.05
mydat.create_artificial({'periods' : 200, 'jump_start' : 46, 'jump_length': 200,'jump_magnitude' : 0.0*noise_magnitude, 'trend_magnitude' : 0.1*noise_magnitude, 'noise' : noise_magnitude})
mydat.resample('Y')
mydat.detect_breakpoints()
print("Now removing the trend through subtracting the linear trend.")
mydat.do_trends()
mydat.remove_trend('linear')
mydat.detect_breakpoints()


# # Example B3:  Breakpoint detection on difference values ("variance")

# In[5]:


mydat = TrendLims1D('De Bilt T2m')
mydat.load_file('../test_data/ehdb_t2m.nc')
# Detect breakpoints
mydat.subset(slice('1950-01-01','2017-12-31'))
mydat.resample('Y')
mydat.remove_trend(method='differences')
mydat.detect_breakpoints()
mydat.plot()


# ## Example B4:  Breakpoint detection on a detrended dataset ("variance")

# In[6]:


mydat = TrendLims1D('De Bilt')
mydat.load_file('../test_data/ehdb_t2m.nc')
mydat.resample('Y')
mydat.subset(slice('1910','2017'))
mydat.do_trends()
mydat.remove_trend('theilsen')
mydat.plot()
mydat.detect_breakpoints()


# ## Example C1: Trend analysis on >100 years of data

# In[7]:


mydat = TrendLims1D('De Bilt')
mydat.load_file('../test_data/ehdb_t2m.nc')
mydat.resample('Y')
mydat.plot()
mydat.do_trends()
mydat.trends


# ## Example C2: Trend analysis on 30 years - Robustness of Theil-Sen estimator explains difference

# In[8]:


mydat.reset()
mydat.resample('Y')
mydat.subset(slice('1987','2018'))
mydat.plot()
mydat.do_trends()
mydat.trends


# ## Example D: Analyse the residuals after subtracting a certain trend

# In[9]:


import matplotlib.pyplot as plt
mydat.reset()
mydat.resample('Y')
mydat.subset(slice('1910','2017'))
mydat.plot('original')
plt.show()
mydat.do_trends()
mydat.remove_trend(method='polynomial')
mydat.plot('detrended')
plt.show()
mydat.do_residual_analysis()
#mydat.do_residual_analysis('linear')
pass # dummy statement for preventing image showing up twice


# In[10]:


mydat.detect_breakpoints()
mydat.breakpoints


# ## Example E: Apply the concept of minimal detectable trends

# In[11]:


mydat.reset()
mydat.weatherhead_framework(trend_magnitude=0.02) # Taken from GCOS table

