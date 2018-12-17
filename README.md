# Trend framework

This is a private repository intended for development of the 'trends and their limits' framework from Oct 2018 within the C3S_511 project.

In the first stage, the diagnostics will be tested on one dimensional data (either station data, or a single gridpoint extracted from a gridded dataset). This allows for the development of robust diagnostics that can later be applied to gridded data. The folder test_data contains a file with 2m temperature from the station 'De Bilt' in the Netherlands. This timeseries is relatively long and is homogenized.

The file _c3s_511_trends.py_ contains the class TrendLims1D in which different diagnostics are implemented, and more will be implemented in the future. The files with the extension `.ipynb` are Jupyter notebook files and serve as examples for using the framework. For users who have no experience with the Jupyter notebook, it is recommended to read sections 1.1 and 3 from [this tutorial](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/what_is_jupyter.html#notebook-document).

# Installation

## Creating the environment
 - Create the environment needed for ESMValTool2 as described [here](https://esmvaltool.readthedocs.io/en/version2_development/user_guide2/index.html#installing-esmvaltool) with the Python 3 environment (note for ETH users: replace `python=3` with `python=3.6`)
 - conda install -c r rpy2
 - conda install jupyter
 - conda install statsmodels
 - install the following R packages from an R terminal opened while the conda environment is activated (recommended):
    - iki.dataclim
      - lubridate
      - zoo
      - climdex.pcic
    - trend
    
 - alternatively from Python while the conda environment is activated (not recommended):
 ```
 import rpy2.robjects.packages as rpackages
 import rpy2.robjects as robjects

 # import R's utility package
 utils = rpackages.importr('utils')

 # select a mirror for R packages
 utils.chooseCRANmirror(ind=1) # select the first mirror in the list

 utils.install_packages('trend')
 if rpackages.isinstalled('trend'):
    print("Succesfully installed the 'trend' R package from CRAN")
 else:
    print("Warning: installation of the 'trend' R package was not succesfull") 
 ```
## Clone this repository
Navigate to the directory where you want the software to be installed and clone this repository. The Jupyter notebook is started with the command `jupyter-notebook`.

# ECV evaluator guidance [in construction]

## Product homogeneity
Run breakpoint tests on the values and the variance of the data product. Is the data marked as suspect? Are the outcomes ain agreement with what you would expect as an expert? Are the findings consistent with the information that is provided by the dataprovider or any other studies that looked into the homogeneity of the product?

## Is the timeseries long enough to assess GCOS requirements?
Apply the weatherhead framework to the data product. Is the dataset long enough to assess the GCOS requirements?

# Usefull material
- ATBD document from ECA&D: https://www.ecad.eu/documents/atbd.pdf
- R-package trend documentation: https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf
- R-package iki.dataclim documentation: https://cran.r-project.org/web/packages/iki.dataclim/iki.dataclim.pdf
- Testing the assumptions of linear regression: http://people.duke.edu/~rnau/testing.htm
- Understanding Q-Q Plots: https://data.library.virginia.edu/understanding-q-q-plots/
- Mann-Kendall test  https://vsp.pnnl.gov/help/vsample/design_trend_mann_kendall.htm
