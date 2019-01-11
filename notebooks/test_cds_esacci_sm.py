
# coding: utf-8

# In[2]:


import sys
sys.path.append("..")
from c3s_511_trends import TrendLims1D


# # Dataset description
# ESA CCI Soilmoisture Combined. A soilmoisture product that blends different satellite products (both active and passive). 
# 
# The selected gridpoint SAUD (22.1N;50.8E) corresponds to desert in Saudi Arabia. Su et al. 2006 [1] shows detected breakpoints in their Figure 2 using two different tests. For the selected gridpoint they find a breakpoint at the start of 1998 (t3,4) and another one halfway 2002 (t4,5) using alpha=0.01 for both tests. As can be seen below, the dataset contains several periods with missing data. We get rid of these missing data through resampling to monthly means as is done in Su et al. 2006. 
# 
# The selected gridpoint AUST (19.1S;132.6E) is in the central northern territory in Australia. Also for this gridpoint, a breakpoint is found by Su et al. 2006.
# 
# [1] Su, C.-H., D. Ryu, W. Dorigo,
# S. Zwieback, A. Gruber, C. Albergel,
# R. H. Reichle, and W. Wagner
# (2016), Homogeneity of a global
# multisatellite soil moisture cli-
# mate data record, Geophys.
# Res. Lett., 43, 11,245–11,252,
# doi:10.1002/2016GL070458.

# # Breakpoint detection
# 

# In[4]:


saud = TrendLims1D('CDS ESACCI SM - SAUD')
saud.load_file('../test_data/cds_esacci_sm_SAUD.nc',var_name='sm')
saud.subset(slice('1988','2017'))
saud.plot('SAUD daily')
saud.resample('M')
saud.plot('SAUD monthly')
saud.breakpoint_recipe_differences(resamplefreq='M')
saud.breakpoint_recipe_values(resamplefreq='M')


# ### The timeseries is indeed marked as suspect for the breakpoint test applied to values and to differences. This is consistent with Su et al. 2006. 

# In[5]:


oman = TrendLims1D('CDS ESACCI SM - OMAN')
oman.load_file('../test_data/cds_esacci_sm_OMAN.nc',var_name='sm')
oman.subset(slice('1988','1997'))
oman.plot('OMAN daily')
oman.resample('M')
oman.plot('OMAN monthly')
saud.breakpoint_recipe_differences(resamplefreq='M')
saud.breakpoint_recipe_values(resamplefreq='M')


# ### The timeseries for OMAN is also marked as suspect for the breakpoint test applied to values and to differences. This is consistent with Su et al. 2006. 

# In[6]:


aust = TrendLims1D('CDS ESACCI SM - AUST')
aust.load_file('../test_data/cds_esacci_sm_AUST.nc',var_name='sm')
aust.subset(slice('1991','2017'))
aust.plot('AUST daily')
aust.resample('M')
aust.plot('AUST monthly')
# Note that now the breakpoint detection is applied to yearly values (the default for these functions), 
# this choice was made since there were data gaps
aust.breakpoint_recipe_differences()
aust.breakpoint_recipe_values()


# In[7]:


# Now try another subset, and use monthly values
aust.reset()
aust.subset(slice('2002','2017'))
aust.plot('AUST daily')
aust.resample('M')
aust.plot('AUST monthly')
# Note that now the breakpoint detection is applied to yearly values (the default for these functions), 
# this choice was made since there were data gaps
aust.breakpoint_recipe_differences(resamplefreq='M')
aust.breakpoint_recipe_values(resamplefreq='M')


# # Timeseries length needed to assess GCOS stability for three different locations

# In[8]:


oman.weatherhead_framework(trend_magnitude=0.001)


# In[9]:


saud.weatherhead_framework(trend_magnitude=0.001)


# In[10]:


aust.weatherhead_framework(trend_magnitude=0.001)


# ### The timeseries length needed to assess GCOS stability differs from ~40 years for SAUD/OMAN to ~97 years for Australia

# # Conclusions
# 
#  - Testing three different grid points in a 'blended' satellite product covering about 30 years
# 
#  - It was found that breakpoints are detected when they are expected to occur according to Su et al. It was not yet tested how the tests perform in regions where Su et al. do not detect any breakpoints
#  
#  - The timeseries length needed to assess GCOS stability shows a strong variation between the different timeseries, from 40 to 97 years. These values seem reasonable.
