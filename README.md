# trend_framework

This is a private repository intended for development of the 'trends and their limits' framework from Oct 2018 within the C3S_511 project.

In the first stage, the diagnostics will be tested on station data. This allows for the development of robust diagnostics that can later be applied to gridded data. The folder test_data contains a file with 2m temperature from the station 'De Bilt' in the Netherlands. This timeseries is relatively long, and is non-homogenized, should therefore contain breakpoints.

The file _c3s_511_trends.py_ contains the class TrendLims1D in which different diagnostics are implemented, and more will be implemented in the future. The file _framework.ipynb_ (iPython Notebook) contains several examples and can be viewed on Github.



# Usefull material
- Testing the assumptions of linear regression: http://people.duke.edu/~rnau/testing.htm
- Understanding Q-Q Plots: https://data.library.virginia.edu/understanding-q-q-plots/
- Mann-Kendall test  https://vsp.pnnl.gov/help/vsample/design_trend_mann_kendall.htm
