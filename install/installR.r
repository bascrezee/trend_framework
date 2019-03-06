# edit the line below to the right directories
# the right directory can be found by activating the conda environment and typing 'which R'
# which returns:  '/net/exo/landclim/crezees/conda/envs/trends/bin/R' in this case
installdir <- "/net/exo/landclim/crezees/conda/envs/trends/lib/R/library/"
install.packages("forecast",installdir,repos='http://cran.us.r-project.org')
install.packages("bfast",installdir,repos='http://cran.us.r-project.org')
install.packages("iki.dataclim",installdir,repos='http://cran.us.r-project.org')
install.packages("trend",installdir,repos='http://cran.us.r-project.org')
