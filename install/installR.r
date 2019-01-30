#!/net/exo/landclim/crezees/conda/envs/trends/bin/Rscript
installdir <- "/net/exo/landclim/crezees/conda/envs/trends/lib/R/library/"
install.packages("forecast",installdir,repos='http://cran.us.r-project.org')
install.packages("bfast",installdir,repos='http://cran.us.r-project.org')
install.packages("iki.dataclim",installdir,repos='http://cran.us.r-project.org')
install.packages("trend",installdir,repos='http://cran.us.r-project.org')
