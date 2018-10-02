# trend_framework

This is a private repository intended for development of the 'trends and their limits' framework from Oct 2018 within the C3S_511 project.

In the first stage, the diagnostics will be tested on station data. This allows for the development of robust diagnostics that can later be applied to gridded data. The folder test_data contains a file with 2m temperature from the station 'De Bilt' in the Netherlands. This timeseries is relatively long, and is non-homogenized, should therefore contain breakpoints.

The file _c3s_511_trends.py_ contains the class TrendLims1D in which different diagnostics are be implemented. iPython notebooks can be used to create different trend analyses, inspect and plot the results.
