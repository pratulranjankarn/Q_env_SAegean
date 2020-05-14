#!/bin/bash

rm -rf Intgrl_files

mkdir Intgrl_files

while read A1 A2 A3;
do
grep "${A1}" intgrl_cycnt_compil_coda_nr.txt > ./Intgrl_files/${A1}_intgrl.txt

done < stlist_cycnt.txt

mkdir station_files

for i in ./Intgrl_files/*
do
j=${i//Intgrl_files}
k=${j//_intgrl.txt}
mw=`grep "4-8 UC_norm" ${i}`
qw=`grep "2-4 UC_norm" ${i}`
rw=`grep "8-16 UC_norm" ${i}`
temp=`grep "${k#.//}" intgrl_egel_compil_coda_nr.txt | tail -n1`
temp2=`grep "${k#.//}" intgrl_husn_compil_coda_nr.txt | tail -n1`
if [[ ${#mw} -gt 0 ]] || [[ ${#qw} -gt 0 ]] || [[ ${#rw} -gt 0 ]]; then
echo ${#mw} ${#qw} ${#rw} ${k#.//}
#: <<EOF
mkdir ./station_files/${k#.//}
#echo $temp
if [[ "${#temp}" -eq 0 ]] && [[ "${#temp2}" -eq 0 ]]
then
echo ${k#.//}
rm -f ./station_files/${k#.//}/U*coda*.txt
fi
grep "1-2 UC_norm" ${i} >> ./station_files/${k#.//}/U1_coda_1-2.txt
grep "2-4 UC_norm" ${i} >> ./station_files/${k#.//}/U1_coda_2-4.txt
grep "4-8 UC_norm" ${i} >> ./station_files/${k#.//}/U1_coda_4-8.txt
grep "8-16 UC_norm" ${i} >> ./station_files/${k#.//}/U1_coda_8-16.txt
grep "16-32 UC_norm" ${i} >> ./station_files/${k#.//}/U1_coda_16-32.txt

grep "1-2 UC_norm" ${i} >> ./station_files/${k#.//}/U2_coda_1-2.txt
grep "2-4 UC_norm" ${i} >> ./station_files/${k#.//}/U2_coda_2-4.txt
grep "4-8 UC_norm" ${i} >> ./station_files/${k#.//}/U2_coda_4-8.txt
grep "8-16 UC_norm" ${i} >> ./station_files/${k#.//}/U2_coda_8-16.txt
grep "16-32 UC_norm" ${i} >> ./station_files/${k#.//}/U2_coda_16-32.txt

grep "1-2 UC_norm" ${i} >> ./station_files/${k#.//}/U3_coda_1-2.txt
grep "2-4 UC_norm" ${i} >> ./station_files/${k#.//}/U3_coda_2-4.txt
grep "4-8 UC_norm" ${i} >> ./station_files/${k#.//}/U3_coda_4-8.txt
grep "8-16 UC_norm" ${i} >> ./station_files/${k#.//}/U3_coda_8-16.txt
grep "16-32 UC_norm" ${i} >> ./station_files/${k#.//}/U3_coda_16-32.txt
#EOF
fi
done

rm -rf Intgrl_files
