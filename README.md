Trend framework
===============

This is a private repository intended for development of the 'trends and their limits' framework from Oct 2018 within the C3S_511 project.

In the first stage, the diagnostics will be tested on one dimensional data (either station data, or a single gridpoint extracted from a gridded dataset). This allows for the development of robust diagnostics that can later be applied to gridded data. The folder test_data contains a file with 2m temperature from the station 'De Bilt' in the Netherlands. This timeseries is relatively long and is homogenized.

The file _c3s_511_trends.py_ contains the class TrendLims1D in which different diagnostics are implemented, and more will be implemented in the future. The files with the extension `.ipynb` are Jupyter notebook files and serve as examples for using the framework. For users who have no experience with the Jupyter notebook, it is recommended to read sections 1.1 and 3 from [this tutorial](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/what_is_jupyter.html#notebook-document).

Installation
===============

### Creating the environment
 - Create the environment needed for ESMValTool2 as described [here](https://esmvaltool.readthedocs.io/en/version2_development/user_guide2/index.html#installing-esmvaltool) with the Python 3 environment (note for ETH users: replace `python=3` with `python=3.6`)
 - conda install -c r rpy2
 - conda install jupyter
 - conda install statsmodels
 - conda install -c conda-forge pathos
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
### Clone this repository
Navigate to the directory where you want the software to be installed and clone this repository. The Jupyter notebook is started with the command `jupyter-notebook`.

Framework testing guidance [in construction]
===========================================
The aim of this framework is to (1) analyze the dataset on the presence of potential breakpoints and (2) to calculate the length of the timeseries that is needed to detect a trend larger than a certain magnitude (which can be taken to be the GCOS stability). The framework is in development, therefore, don't hesitate to provide feedback, contribute to the code, or to the documentation by opening an issue.

### Download test data for your ECV from the CDS
Use the script 'cds_toolbox_script.py' as an example script for downloading data from the Climate Data Store (CDS). It can be copy-pasted into the Toolbox and adapted towards your own needs. Choose 2-3 different grid points to do the analysis.

### Perform the analysis
Insert the Python code from the file 'test_cds_esacci_sm.py' into a jupyter notebook. Follow the given structure, thereby answering the following questions:

##### Dataset description
Give a brief description of the dataset focussing on information that might be relevant for the trend framework. Briefly describe which locations you selected to do the analysis.

##### Product homogeneity
Run breakpoint tests on the values and the variance of the data product. Is the data marked as suspect? Are the outcomes in agreement with what you would expect as an expert? Are the findings consistent with the information that is provided by the dataprovider or any other studies that looked into the homogeneity of the product?

##### Is the timeseries long enough to assess the GCOS stability requirements?
Apply the weatherhead framework to the data product. Is the dataset long enough to assess the GCOS requirements?

### Upload your test results to this repository to share with others
Once you have completed the testing, export your notebook file as a PDF. Upload the PDF to the repository in the directory PDFs.
First make sure your local repository is up-to-date with the remote repository by executing
```
git pull
```
See an overview of the files that you have added yourself
```
git status
```
Add all the files that you have added
```
git add myfile
```
Commit your results (this is stil a local command) with a usefull commit message
```
git commit -m "Added PDF test report for soil moisture"
```
Push your results to the remote repository
```
git push
```

Trends and their limits evaluator guidance [in construction]
============================================================

Detecting trends and breakpoints in timeseries is a challenging task. Not the least, since the two phenomena are linked, and therefore a priori knowledge is needed for a correct approach. For this reason, the evaluator guidance is split into two parts, where the first part is an exploratory phase, where you will evaluate in-depth two or three gridpoints from your dataset. 

#### Exploratory phase
Select two or three gridpoints from your dataset. Think about a good temporal coverage, and in your selection try to take gridpoints that are representative for a larger subset of the dataset. Also, the selection can take into account earlier studies, that allow for comparison of the results.

##### Trends
Apply the trend tests to the data. Is there a trend present in the data? Do the Theil-Sen slope estimator and the linear regression estimate of the trend agree? Is the trend constant over time? Or can you identify subsets in the data for which the trend is different. If so, briefly describe the nature of the trend.

##### Breakpoints
If there is a trend present in the data, it needs to be removed before breakpoint detection is applied. Based on your prior knowledge of the trend within the data, remove it either using the Theil-Sen slope estimate (trend is constant over time) or the polynomial fit (trend changes over time).

#### Creation of maps [in construction]
 - create your own trend and lims recipe (using knowledge from exploratory phase) and run it on the gridded dataset at yearly resolution
 - produce maps with:
     - Theil-Sen slope estimate
     - Mann-Kendall p-value categories (p<0.05 ; p<0.01)
     - breakpoint category (data is: usefull, doubtfull, suspect)
     - timestamps of breakpoints that are abundant
      
     




# Usefull material
- ATBD document from ECA&D: https://www.ecad.eu/documents/atbd.pdf
- R-package trend documentation: https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf
- R-package iki.dataclim documentation: https://cran.r-project.org/web/packages/iki.dataclim/iki.dataclim.pdf
- Testing the assumptions of linear regression: http://people.duke.edu/~rnau/testing.htm
- Understanding Q-Q Plots: https://data.library.virginia.edu/understanding-q-q-plots/
- Mann-Kendall test  https://vsp.pnnl.gov/help/vsample/design_trend_mann_kendall.htm
