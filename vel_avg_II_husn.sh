#!/bin/bash

for i in ./station_files/*
do
echo ${i:16:4}
aw=`find ./husn_events/ -name "husn.20*.pick" -exec grep "${i:16:4} " {} + | grep "S_"`
if [[ "${#aw}" -gt 0 ]]; then
find ./husn_events/ -name "husn.20*.pick" -exec grep "${i:16:4} " {} + | grep "S_" > temp.txt
#rm -f ${i}/vel.txt ${i}/vel_avg_husn.txt
while read V1 V2 V3 V4 V5 V6 V7 V8; 
do
Otmstmp=`grep -A1 "Date" ${V1:0:59} | tail -n1 | awk '{ print $2 }'`
Odep=`grep -A1 "Date" ${V1:0:59} | tail -n1 | awk '{print $9}'`
Oyeer=`grep -A1 "Date" ${V1:0:59} | tail -n1 | awk '{ print $1 }'`
if [[ "${Oyeer:0:4}" -lt 2019 ]]; then
Otmstmp=${Otmstmp}"0"
#echo $Otmstmp
Stmstmp=${V5}
#echo $Stmstmp
epdist=${V2}
distnce=$(echo "sqrt(${Odep}^2 + ${epdist}^2)" | bc -l)
if (( $(echo "${distnce} <= 100" | bc -l) )); then
echo ${V1:0:59}
if [[ "${Otmstmp:3:2}" == "${Stmstmp:3:2}" ]]; then
trvltm=$(echo "${Stmstmp:6:6} - ${Otmstmp:6:6}" | bc -l)
else
trvltm=$(echo "${Stmstmp:6:6} + 60 - ${Otmstmp:6:6}" | bc -l)
fi
echo $Otmstmp $Stmstmp $trvltm
echo ${distnce} success
echo "${distnce}*1000/${trvltm}" | bc -l
echo "${distnce}*1000/${trvltm}" | bc -l >> ${i}/vel.txt
fi
fi
done < temp.txt
#: <<KOD
python <<EOF
import numpy as np
a = np.loadtxt('${i}/vel.txt')
f = open('${i}/vel_avg_husn.txt','w')
f.write("%f" % np.mean(a))
f.close()
quit()
EOF
#KOD
fi
done
