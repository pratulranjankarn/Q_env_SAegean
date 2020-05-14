#!/bin/bash

for i in ./station_files/*
do
echo ${i:16:4}
aw=`find ./husn_events/ -name "husn.201101*.pick" -exec grep " ${i:16:4} " {} +`
if [[ "${#aw}" -gt 0 ]]; then
find ./husn_events/ -name "husn.201101*.pick" -exec grep " ${i:16:4} " {} + > temp.txt
#rm -f ${i}/vel.txt ${i}/vel_avg_husn.txt
while read V1 V2 V3 V4 V5 V6 V7 V8; 
do 
#echo $V1;
if [[ ${#V4} -lt 10 ]]; then
Ptmstmp=${V4}"0"${V5}
distnce=${V6}
else
Ptmstmp=${V4}
distnce=${V5}
fi
if (( $(echo "${distnce} <= 100" | bc -l) )); then
Stmstmp=`grep -A1 " ${V2}" ${V1%:} | tail -n1 | grep " ES" | awk '{ print $2,$3,$4}' | awk '{
if (length($3) >= 2)
print $1"0"$2;
else
print $1;
}';`
Otmstmp=`grep " H=" ${V1%:} | awk '{ print $2 }'`
if [[ "${Otmstmp: -4:1}" == ":" ]]; then
Otmstmp=${Otmstmp:2:6}"0"${Otmstmp: -3:3}
else
Otmstmp=${Otmstmp:2:6}${Otmstmp: -4:4}
fi
#echo ${Otmstmp}
if [[ "${Otmstmp:3:2}" == "${Ptmstmp:3:2}" ]]; then
trvltm=$(echo "${Ptmstmp:6:4} - ${Otmstmp:6:4}" | bc -l)
else
trvltm=$(echo "${Ptmstmp:6:4} + 60 - ${Otmstmp:6:4}" | bc -l)
fi
#echo $trvltm
if [[ "${#Stmstmp}" -gt 0 ]]; then
echo $V1
echo ${Otmstmp}
if [[ "${#Stmstmp}" -eq 5 ]]; then
#echo ${Stmstmp}
timediff=$(echo "${Stmstmp: -4:4} - ${Ptmstmp: -4:4}" | bc -l)
elif [[ "${#Stmstmp}" -eq 6 ]]; then
timediff=$(echo "${Stmstmp: -3:3} + 60 - ${Ptmstmp: -4:4}" | bc -l)
else
timediff=$(echo "${Stmstmp: -4:4} + 60 - ${Ptmstmp: -4:4}" | bc -l)
fi
Strvltm=$(echo "${trvltm} + ${timediff}" | bc -l)
echo $Ptmstmp $Stmstmp $timediff $Strvltm

echo ${distnce} success
echo "${distnce}*1000/${Strvltm}" | bc -l
echo "${distnce}*1000/${Strvltm}" | bc -l >> ${i}/vel.txt
fi
fi
done < temp.txt
python <<EOF
import numpy as np
a = np.loadtxt('${i}/vel.txt')
f = open('${i}/vel_avg_husn.txt','w')
f.write("%f" % np.mean(a))
f.close()
quit()
EOF
fi
done
