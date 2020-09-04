#!/bin/bash


branchname=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')

trendcli -p '/net/exo/landclim/data/dataset/ERA5_deterministic/recent/0.25deg_lat-lon_1m/processed/regrid/era5_deterministic_recent.swvl1.025deg.1m.2???.nc' -v swvl1 --datasetname branch-${branchname}_ERA5
