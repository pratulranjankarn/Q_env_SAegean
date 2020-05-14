#!/bin/bash

rm -f comb_12giwr_2s_100_intgrl.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_1-2.txt" -size +0c -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_1-2.txt" -size +0c -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_1-2.txt" -size +0c -exec grep " " {} + >> check.txt

while read stnm A2 A3;
do
#stnm="ANAF"
echo ${stnm}

grep -w "${stnm}" check.txt > check2.txt 

if [ -a ./station_files/${stnm}/vel_avg_husn.txt ]; then
v_avg=`cat ./station_files/${stnm}/vel_avg_husn.txt`
else
v_avg="3300"
fi
rm -f  ${stnm}_g.txt

IFS=' '

while read V1 V2 V3 V4 V5 V6;
do
IFS='/'
read -ra arr <<< "${V1}"
IFS=' '
orign=${arr[2]}
#echo $orign
if [ "${orign:0:1}" = "0" ]; then
yeer="20"${orign:0:2}
mnth=${orign:2:2}
else
yeer=${orign:0:4}
mnth=${orign:5:2}
fi
temp=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$V3"`
#echo ${arr[0]} ${arr[1]} ${arr[2]}
if (( $(echo "${V2} < 100" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then

sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "*.hyp"`
read -ra hyp_ar <<< "${sw}"
if [ "${#hyp_ar[1]}" -gt 0 ]; then
echo ${#hyp_ar[1]}
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_12giwr_2s_100_intgrl.txt
IFS=' '
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_12giwr_2s_100_intgrl.txt
IFS=' '
fi
fi
done < check2.txt

done < stlist_f.txt




rm -f comb_24giwr_2s_100_intgrl.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_2-4.txt" -size +0c -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_2-4.txt" -size +0c -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_2-4.txt" -size +0c -exec grep " " {} + >> check.txt

while read stnm A2 A3;
do
#stnm="ANAF"
echo ${stnm}

grep -w "${stnm}" check.txt > check2.txt 

if [ -a ./station_files/${stnm}/vel_avg_husn.txt ]; then
v_avg=`cat ./station_files/${stnm}/vel_avg_husn.txt`
else
v_avg="3300"
fi
rm -f  ${stnm}_g.txt

IFS=' '

while read V1 V2 V3 V4 V5 V6;
do
IFS='/'
read -ra arr <<< "${V1}"
IFS=' '
orign=${arr[2]}
#echo $orign
if [ "${orign:0:1}" = "0" ]; then
yeer="20"${orign:0:2}
mnth=${orign:2:2}
else
yeer=${orign:0:4}
mnth=${orign:5:2}
fi
temp=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$V3"`
#echo ${arr[0]} ${arr[1]} ${arr[2]}
if (( $(echo "${V2} < 100" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then

sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "*.hyp"`
read -ra hyp_ar <<< "${sw}"
if [ "${#hyp_ar[1]}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_24giwr_2s_100_intgrl.txt
IFS=' '
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_24giwr_2s_100_intgrl.txt
IFS=' '
fi
fi
done < check2.txt

done < stlist_f.txt




rm -f comb_48giwr_2s_100_intgrl.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_4-8.txt" -size +0c -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_4-8.txt" -size +0c -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_4-8.txt" -size +0c -exec grep " " {} + >> check.txt

while read stnm A2 A3;
do
#stnm="ANAF"
echo ${stnm}

grep -w "${stnm}" check.txt > check2.txt 

if [ -a ./station_files/${stnm}/vel_avg_husn.txt ]; then
v_avg=`cat ./station_files/${stnm}/vel_avg_husn.txt`
else
v_avg="3300"
fi
rm -f  ${stnm}_g.txt

IFS=' '

while read V1 V2 V3 V4 V5 V6;
do
IFS='/'
read -ra arr <<< "${V1}"
IFS=' '
orign=${arr[2]}
#echo $orign
if [ "${orign:0:1}" = "0" ]; then
yeer="20"${orign:0:2}
mnth=${orign:2:2}
else
yeer=${orign:0:4}
mnth=${orign:5:2}
fi
temp=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$V3"`
#echo ${arr[0]} ${arr[1]} ${arr[2]}
if (( $(echo "${V2} < 100" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then

sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "*.hyp"`
read -ra hyp_ar <<< "${sw}"
if [ "${#hyp_ar[1]}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_48giwr_2s_100_intgrl.txt
IFS=' '
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_48giwr_2s_100_intgrl.txt
IFS=' '
fi
fi
done < check2.txt

done < stlist_f.txt



rm -f comb_816giwr_2s_100_intgrl.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_8-16.txt" -size +0c -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_8-16.txt" -size +0c -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_8-16.txt" -size +0c -exec grep " " {} + >> check.txt

while read stnm A2 A3;
do
#stnm="ANAF"
echo ${stnm}

grep -w "${stnm}" check.txt > check2.txt 

if [ -a ./station_files/${stnm}/vel_avg_husn.txt ]; then
v_avg=`cat ./station_files/${stnm}/vel_avg_husn.txt`
else
v_avg="3300"
fi
rm -f  ${stnm}_g.txt

IFS=' '

while read V1 V2 V3 V4 V5 V6;
do
IFS='/'
read -ra arr <<< "${V1}"
IFS=' '
orign=${arr[2]}
#echo $orign
if [ "${orign:0:1}" = "0" ]; then
yeer="20"${orign:0:2}
mnth=${orign:2:2}
else
yeer=${orign:0:4}
mnth=${orign:5:2}
fi
temp=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$V3"`
#echo ${arr[0]} ${arr[1]} ${arr[2]}
if (( $(echo "${V2} < 100" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then

sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "*.hyp"`
read -ra hyp_ar <<< "${sw}"
if [ "${#hyp_ar[1]}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_816giwr_2s_100_intgrl.txt
IFS=' '
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_816giwr_2s_100_intgrl.txt
IFS=' '
fi
fi
done < check2.txt

done < stlist_f.txt



rm -f comb_1632giwr_2s_100_intgrl.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_16-32.txt" -size +0c -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_16-32.txt" -size +0c -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_16-32.txt" -size +0c -exec grep " " {} + >> check.txt

while read stnm A2 A3;
do
#stnm="ANAF"
echo ${stnm}

grep -w "${stnm}" check.txt > check2.txt 

if [ -a ./station_files/${stnm}/vel_avg_husn.txt ]; then
v_avg=`cat ./station_files/${stnm}/vel_avg_husn.txt`
else
v_avg="3300"
fi
rm -f  ${stnm}_g.txt

IFS=' '

while read V1 V2 V3 V4 V5 V6;
do
IFS='/'
read -ra arr <<< "${V1}"
IFS=' '
orign=${arr[2]}
#echo $orign
if [ "${orign:0:1}" = "0" ]; then
yeer="20"${orign:0:2}
mnth=${orign:2:2}
else
yeer=${orign:0:4}
mnth=${orign:5:2}
fi
temp=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$V3"`
#echo ${arr[0]} ${arr[1]} ${arr[2]}
if (( $(echo "${V2} < 100" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then

sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "*.hyp"`
read -ra hyp_ar <<< "${sw}"
if [ "${#hyp_ar[1]}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/cyc*.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_1632giwr_2s_100_intgrl.txt
IFS=' '
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_1632giwr_2s_100_intgrl.txt
IFS=' '
fi
fi
done < check2.txt

done < stlist_f.txt

