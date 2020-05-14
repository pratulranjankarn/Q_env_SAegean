#!/bin/bash

for i in /home/user/env_EGEL/husn_events/{2013..2019..1}*
do
cd $i
echo $i
find ./ -name "*..HHZ.D.SAC" > /home/user/env_EGEL/trytst.txt
while read V1
do
alpha=${V1//..HHZ.D.SAC}
beta=${alpha: -4:4}
gamma=${beta//.}
echo $gamma
#: <<SOD
aw=`grep "${gamma}" ./env_integral.txt | grep "1-2" | grep "U1_norm" | awk '{ print $3 }'`
if [[ "${#aw}" -eq 0 ]]; then
rm -f *${gamma}*1-2.sacii
fi
bw=`grep "${gamma}" ./env_integral.txt | grep "2-4" | grep "U1_norm" | awk '{ print $3 }'`
if [[ "${#bw}" -eq 0 ]]; then
rm -f *${gamma}*2-4.sacii
fi
cw=`grep "${gamma}" ./env_integral.txt | grep "4-8" | grep "U1_norm" | awk '{ print $3 }'`
if [[ "${#cw}" -eq 0 ]]; then
rm -f *${gamma}*4-8.sacii
fi
dw=`grep "${gamma}" ./env_integral.txt | grep "8-16" | grep "U1_norm" | awk '{ print $3 }'`
if [[ "${#dw}" -eq 0 ]]; then
rm -f *${gamma}*8-16.sacii
fi
ew=`grep "${gamma}" ./env_integral.txt | grep "16-32" | grep "U1_norm" | awk '{ print $3 }'`
if [[ "${#ew}" -eq 0 ]]; then
rm -f *${gamma}*16-32.sacii
fi
if [[ "${#aw}" -eq 0 ]] && [[ "${#bw}" -eq 0 ]] && [[ "${#cw}" -eq 0 ]] && [[ "${#dw}" -eq 0 ]] && [[ "${#ew}" -eq 0 ]]; then
rm -f *${gamma}..HHZ.D.SAC *${gamma}..HHE.D.SAC *${gamma}..HHN.D.SAC
fi
#SOD
done < /home/user/env_EGEL/trytst.txt
#rm -f *HH*.sacii *HH*.ascii
cd /home/user/env_EGEL/husn_events/
pwd
done
