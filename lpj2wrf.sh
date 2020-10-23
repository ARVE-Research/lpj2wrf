#!/bin/bash

cycle=cycle06

#------
# 54km

datapath=/soge-home/users/cenv0668/arve/science_projects/WRF_LGM/$cycle/54km

# generate mean temperature

ncks -O --no-abc -d time,0870-01-01,1019-12-31 -v tmp $datapath/LGMclimate.nc $datapath/tmp_last150.nc
cdo timmean $datapath/tmp_last150.nc $datapath/tann_mean_last150.nc

#--

xlen=83
ylen=84

outfile=WRF-LPJ_europe_54km_$cycle.nc

echo $xlen $ylen

sed -e "s/XLEN/$xlen/g" -e "s/YLEN/$ylen/g" postprocess.cdl | ncgen -k nc4 -o $outfile

sed 's/CYCLE/'$cycle'/g' europe_54km.namelist > tmp.namelist

./postprocess tmp.namelist $outfile

rm tmp.namelist

#------
# 18km

datapath=/soge-home/users/cenv0668/arve/science_projects/WRF_LGM/$cycle/18km

# generate mean temperature

ncks -O --no-abc -d time,0870-01-01,1019-12-31 -v tmp $datapath/LGMclimate.nc $datapath/tmp_last150.nc
cdo timmean $datapath/tmp_last150.nc $datapath/tann_mean_last150.nc

#--

xlen=186
ylen=147

outfile=WRF-LPJ_europe_18km_$cycle.nc

echo $xlen $ylen

sed -e "s/XLEN/$xlen/g" -e "s/YLEN/$ylen/g" postprocess.cdl | ncgen -k nc4 -o $outfile

sed 's/CYCLE/'$cycle'/g' europe_18km.namelist > tmp.namelist

./postprocess tmp.namelist $outfile

rm tmp.namelist
