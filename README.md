Trend framework
===============

ECV evaluator guidance
===============================
This is a repository containing some trend functionalities for use within the C3S_511 project. The main functionalities are found in the file `trends_core3d.py` and their usage is shown in the example notebook `example_trends_3d_cds-satellite-soil-moisture.ipynb` that can be found in the directory `notebook`. For users who have no experience with the Jupyter notebook, it is recommended to read sections 1.1 and 3 from [this tutorial](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/what_is_jupyter.html#notebook-document).

Installation for ECV evaluators
===============================
 - Navigate to the directory where you want the software to be installed and clone this repository.
 - Navigate to the `install` directory
 - Create a new conda environment named `trend` (as specified in environment.yml): `conda env create -f environment.yml`
 - Activate the environment: `conda activate trend`
 - Open the example notebook and modify the evaluator specific part for reading in and preprocessing your data

Installation for developers
===========================

 - Navigate to the directory where you want the software to be installed and clone this repository.
 - Navigate to the `install` directory
 - Create a new conda environment named `trend` (as specified in environment.yml): `conda env create -f environment.yml`
 - Activate the environment: `conda activate trend`
 - Add rpy2 as package: `conda install rpy2`
 - Edit the file `installR.r` to point to the right installation path (see file) 
 - Run `Rscript installR.r`
 - Wait for the script to be finished
 - Check the installation by running  the example notebook.

Contributing
=============
![Gitflow Workflow (Copyright: Atlassian)](figs/gitflow.png)

# Usefull material
- ATBD document from ECA&D: https://www.ecad.eu/documents/atbd.pdf
- R-package trend documentation: https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf
- R-package iki.dataclim documentation: https://cran.r-project.org/web/packages/iki.dataclim/iki.dataclim.pdf
- Testing the assumptions of linear regression: http://people.duke.edu/~rnau/testing.htm
- Understanding Q-Q Plots: https://data.library.virginia.edu/understanding-q-q-plots/
- Mann-Kendall test  https://vsp.pnnl.gov/help/vsample/design_trend_mann_kendall.htm
