#!/bin/bash

rm -rf Intgrl_files

mkdir Intgrl_files

stnm="PE01"

#while read A1 A2 A3;
#do
grep "${stnm}" intgrl_egel_Stest.txt > ./Intgrl_files/${stnm}_intgrl.txt

#done < stlist_cycnt.txt

mkdir station_files

for i in ./Intgrl_files/*
do
j=${i//Intgrl_files}
k=${j//_intgrl.txt}
mw=`grep "4-8 U1_norm" ${i}`
qw=`grep "2-4 U1_norm" ${i}`
rw=`grep "8-16 U1_norm" ${i}`
#temp=`grep "${k#.//}" intgrl_egel_compil_sc_nr.txt | tail -n1`
#temp2=`grep "${k#.//}" intgrl_husn_compil_sc_nr.txt | tail -n1`
if [[ ${#mw} -gt 0 ]] || [[ ${#qw} -gt 0 ]] || [[ ${#rw} -gt 0 ]]; then
echo ${#mw} ${#qw} ${#rw} ${k#.//}
#: <<EOF
mkdir ./station_files/${k#.//}
#echo $temp
#if [[ "${#temp}" -eq 0 ]] && [[ "${#temp2}" -eq 0 ]]; then
echo ${k#.//}
rm -f ./station_files/${k#.//}/U*.txt
#fi
grep "1-2 U1_norm" ${i} >> ./station_files/${k#.//}/U1_1-2.txt
grep "2-4 U1_norm" ${i} >> ./station_files/${k#.//}/U1_2-4.txt
grep "4-8 U1_norm" ${i} >> ./station_files/${k#.//}/U1_4-8.txt
grep "8-16 U1_norm" ${i} >> ./station_files/${k#.//}/U1_8-16.txt
grep "16-32 U1_norm" ${i} >> ./station_files/${k#.//}/U1_16-32.txt

grep "1-2 U2_norm" ${i} >> ./station_files/${k#.//}/U2_1-2.txt
grep "2-4 U2_norm" ${i} >> ./station_files/${k#.//}/U2_2-4.txt
grep "4-8 U2_norm" ${i} >> ./station_files/${k#.//}/U2_4-8.txt
grep "8-16 U2_norm" ${i} >> ./station_files/${k#.//}/U2_8-16.txt
grep "16-32 U2_norm" ${i} >> ./station_files/${k#.//}/U2_16-32.txt

grep "1-2 U3_norm" ${i} >> ./station_files/${k#.//}/U3_1-2.txt
grep "2-4 U3_norm" ${i} >> ./station_files/${k#.//}/U3_2-4.txt
grep "4-8 U3_norm" ${i} >> ./station_files/${k#.//}/U3_4-8.txt
grep "8-16 U3_norm" ${i} >> ./station_files/${k#.//}/U3_8-16.txt
grep "16-32 U3_norm" ${i} >> ./station_files/${k#.//}/U3_16-32.txt
#EOF
fi
done

rm -rf Intgrl_files
