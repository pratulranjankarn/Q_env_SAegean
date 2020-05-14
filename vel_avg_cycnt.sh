#!/bin/bash

for i in ./station_files/*
do
echo ${i:16:4}
aw=`find ./Cycnet_events/ -name "*.hyp" -exec grep "${i:16:4} " {} + | grep "? S      ?"`
if [[ "${#aw}" -gt 0 ]]; then
find ./Cycnet_events/ -name "*.hyp" -exec grep "${i:16:4} " {} + | grep "? S      ?" > temp.txt

mw=`grep "${i:16:4}" ./stlist_egel.txt`
qw=`grep "${i:16:4}" ./stlist_husn.txt`
sw=`grep "${i:16:4}" ./stlist_cycnt.txt`

if [[ "${#sw}" -gt 0 ]] && [[ "${#qw}" -eq 0 ]] && [[ "${#mw}" -eq 0 ]]; then
rm -f ${i}/vel.txt ${i}/vel_avg_husn.txt
fi

while read V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21 V22 V23 V24 V25 V26 V27; 
do
Odep=`grep "GEOGRAPHIC" ${V1:0:45} | awk '{ print $14 }'`
epdist=${V22}
distnce=$(echo "sqrt(${Odep}^2 + ${epdist}^2)" | bc -l)
if (( $(echo "${distnce} <= 100" | bc -l) )); then
echo ${V1:0:45}
Ohrmin=`grep "GEOGRAPHIC" ${V1:0:45} | awk '{ print $6$7 }'`
Osec=`grep "GEOGRAPHIC" ${V1:0:45} | awk '{ print $8 }'`
echo $Ohrmin $Osec
Shrmin=${V8}
Ssec=${V9}
echo $Shrmin $Ssec
trvltm=$(echo "${V16} + ${V17}" | bc -l)
echo $trvltm
#: <<MOD
echo ${distnce} success
echo "${distnce}*1000/${trvltm}" | bc -l
echo "${distnce}*1000/${trvltm}" | bc -l >> ${i}/vel.txt
#MOD
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
