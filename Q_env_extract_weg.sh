#!/bin/bash

rm -f comb_12giwr_weg.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_1-2.txt" -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_1-2.txt" -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_1-2.txt" -exec grep " " {} + >> check.txt

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
if (( $(echo "${V2} < 60" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then
if [ \( "${mnth}" -le 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -eq 2010 ]; then
#ls ${arr[0]}/${arr[1]}/${arr[2]}/*.pick
evlat=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $2 }'`
evlat="$(echo -e "${evlat}" | tr -d '[:space:]')"
evlat=${evlat//N} 
evlong=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $4 }'`
evlong="$(echo -e "${evlong}" | tr -d '[:space:]')"
evlong=${evlong//E} 
evdep=`grep -w "LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $6 }'`
evdep="$(echo -e "${evdep}" | tr -d '[:space:]')"
evdep=${evdep//Km} 
elif [ \( "${mnth}" -gt 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -gt 2011 ]; then
evlat=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $4 }'` 
evlong=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $5 }'`
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $6 }'`
if [ "${#evdep}" -ge 3 ]; then
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $9 }'`
fi
#echo ${arr[0]}/${arr[1]}/${arr[2]} $evlat $evlong $evdep
else
sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "0*.hyp"`
if [ "${#sw}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $14 }'`
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
fi
fi
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_12giwr_weg.txt
fi 
IFS=' '
done < check2.txt

done < stlist_f.txt




rm -f comb_24giwr_weg.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_2-4.txt" -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_2-4.txt" -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_2-4.txt" -exec grep " " {} + >> check.txt

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
if (( $(echo "${V2} < 60" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then
if [ \( "${mnth}" -le 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -eq 2010 ]; then
#ls ${arr[0]}/${arr[1]}/${arr[2]}/*.pick
evlat=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $2 }'`
evlat="$(echo -e "${evlat}" | tr -d '[:space:]')"
evlat=${evlat//N} 
evlong=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $4 }'`
evlong="$(echo -e "${evlong}" | tr -d '[:space:]')"
evlong=${evlong//E} 
evdep=`grep -w "LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $6 }'`
evdep="$(echo -e "${evdep}" | tr -d '[:space:]')"
evdep=${evdep//Km} 
elif [ \( "${mnth}" -gt 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -gt 2011 ]; then
evlat=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $4 }'` 
evlong=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $5 }'`
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $6 }'`
if [ "${#evdep}" -ge 3 ]; then
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $9 }'`
fi
#echo ${arr[0]}/${arr[1]}/${arr[2]} $evlat $evlong $evdep
else
sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "0*.hyp"`
if [ "${#sw}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $14 }'`
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
fi
fi
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_24giwr_weg.txt
fi 
IFS=' '
done < check2.txt

done < stlist_f.txt




rm -f comb_48giwr_weg.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_4-8.txt" -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_4-8.txt" -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_4-8.txt" -exec grep " " {} + >> check.txt

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
if (( $(echo "${V2} < 60" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then
if [ \( "${mnth}" -le 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -eq 2010 ]; then
#ls ${arr[0]}/${arr[1]}/${arr[2]}/*.pick
evlat=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $2 }'`
evlat="$(echo -e "${evlat}" | tr -d '[:space:]')"
evlat=${evlat//N} 
evlong=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $4 }'`
evlong="$(echo -e "${evlong}" | tr -d '[:space:]')"
evlong=${evlong//E} 
evdep=`grep -w "LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $6 }'`
evdep="$(echo -e "${evdep}" | tr -d '[:space:]')"
evdep=${evdep//Km} 
elif [ \( "${mnth}" -gt 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -gt 2011 ]; then
evlat=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $4 }'` 
evlong=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $5 }'`
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $6 }'`
if [ "${#evdep}" -ge 3 ]; then
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $9 }'`
fi
#echo ${arr[0]}/${arr[1]}/${arr[2]} $evlat $evlong $evdep
else
sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "0*.hyp"`
if [ "${#sw}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $14 }'`
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
fi
fi
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_48giwr_weg.txt
fi 
IFS=' '
done < check2.txt

done < stlist_f.txt



rm -f comb_816giwr_weg.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_8-16.txt" -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_8-16.txt" -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_8-16.txt" -exec grep " " {} + >> check.txt

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
if (( $(echo "${V2} < 60" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then
if [ \( "${mnth}" -le 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -eq 2010 ]; then
#ls ${arr[0]}/${arr[1]}/${arr[2]}/*.pick
evlat=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $2 }'`
evlat="$(echo -e "${evlat}" | tr -d '[:space:]')"
evlat=${evlat//N} 
evlong=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $4 }'`
evlong="$(echo -e "${evlong}" | tr -d '[:space:]')"
evlong=${evlong//E} 
evdep=`grep -w "LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $6 }'`
evdep="$(echo -e "${evdep}" | tr -d '[:space:]')"
evdep=${evdep//Km} 
elif [ \( "${mnth}" -gt 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -gt 2011 ]; then
evlat=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $4 }'` 
evlong=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $5 }'`
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $6 }'`
if [ "${#evdep}" -ge 3 ]; then
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $9 }'`
fi
#echo ${arr[0]}/${arr[1]}/${arr[2]} $evlat $evlong $evdep
else
sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "0*.hyp"`
if [ "${#sw}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $14 }'`
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
fi
fi
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_816giwr_weg.txt
fi 
IFS=' '
done < check2.txt

done < stlist_f.txt



rm -f comb_1632giwr_weg.txt
rm -f comb_giwr.txt

find ./husn_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_16-32.txt" -exec grep " " {} + > check.txt
find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_16-32.txt" -exec grep " " {} + >> check.txt
find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "env_param_weg_16-32.txt" -exec grep " " {} + >> check.txt

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
if (( $(echo "${V2} < 60" | bc -l) )) && (( $(echo "${temp} < 0.1" | bc -l) )) && (( $(echo "${temp} > 0.00001" | bc -l) )) && (( $(echo "${V4} > 0" | bc -l) )); then
if [ \( "${mnth}" -le 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -eq 2010 ]; then
#ls ${arr[0]}/${arr[1]}/${arr[2]}/*.pick
evlat=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $2 }'`
evlat="$(echo -e "${evlat}" | tr -d '[:space:]')"
evlat=${evlat//N} 
evlong=`grep -w " LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $4 }'`
evlong="$(echo -e "${evlong}" | tr -d '[:space:]')"
evlong=${evlong//E} 
evdep=`grep -w "LAT=" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | awk '{ print $6 }'`
evdep="$(echo -e "${evdep}" | tr -d '[:space:]')"
evdep=${evdep//Km} 
elif [ \( "${mnth}" -gt 01 -a "${yeer}" -eq 2011 \) -o "${yeer}" -gt 2011 ]; then
evlat=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $4 }'` 
evlong=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $5 }'`
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $6 }'`
if [ "${#evdep}" -ge 3 ]; then
evdep=`grep -w -A1 "Latitude Longitude" ${arr[0]}/${arr[1]}/${arr[2]}/*.pick | tail -n1 | awk '{ print $9 }'`
fi
#echo ${arr[0]}/${arr[1]}/${arr[2]} $evlat $evlong $evdep
else
sw=`find ${arr[0]}/${arr[1]}/${arr[2]}/ -name "0*.hyp"`
if [ "${#sw}" -gt 0 ]; then
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/0*.hyp | awk '{ print $14 }'`
else
evlat=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $10 }'` 
evlong=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $12 }'`
evdep=`grep -w "GEOGRAPHIC" ${arr[0]}/${arr[1]}/${arr[2]}/*grid0.loc.hyp | awk '{ print $14 }'`
fi
fi
echo ${arr[2]} ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6}
stltln=`grep -w "${stnm}" stlist_f.txt | awk '{ print $2,$3 }'`
echo ${evlat} ${evlong} ${evdep} ${stltln} ${v_avg} ${V2} ${V3} ${V4} ${V5} ${V6} >> comb_1632giwr_weg.txt
fi 
IFS=' '
done < check2.txt

done < stlist_f.txt
