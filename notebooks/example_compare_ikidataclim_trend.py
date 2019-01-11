
# coding: utf-8

# In[1]:


from c3s_511_trends import TrendLims1D


# In[2]:


mydat = TrendLims1D('test data')
mydat.create_artificial()
mydat.plot()
mydat.breakpoint_recipe_values()
mydat.breakpoints


# In[3]:


mydat.breakpoint_recipe_values(rpackagename='iki.dataclim')
mydat.breakpoints

