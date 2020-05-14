#!/bin/bash

: <<SOD

for i in ~/env_EGEL/station_files/*
do
echo $i
cp ~/env_EGEL/station_files/TUR1/intgrl_plot.gmt ${i}/
done

for i in ~/env_EGEL/station_files/*
do
echo $i
cd $i
./intgrl_plot.gmt
cd ..
done

SOD

awk '{ print $1 }' stlist_egel.txt > Stest.txt
#echo "NAXO" > Stest.txt

while read i
do
echo ~/env_EGEL/station_files/${i}/
cp ~/env_EGEL/station_files/TUR1/g_func.m ~/env_EGEL/station_files/TUR1/optim_param.m ~/env_EGEL/station_files/TUR1/g_func_plot.m ~/env_EGEL/station_files/${i}/
cd ~/env_EGEL/station_files/${i}/
matlab -nosplash -nodesktop -noFigureWindows <<EOF
optim_param
g_func_plot
quit
EOF
cd ..
done < Stest.txt
