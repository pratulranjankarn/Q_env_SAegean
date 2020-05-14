#!/bin/bash

#find ./EGEL_events/ -name "env_integral.txt" -size +1b -exec grep "PE01" {} + > test.txt
#matlab <<EOF
#fid = fopen('test.txt','r');
#A = textscan(fid,'%s %s %s %s %s');
#fclose(fid);
#B = string(A{1});
#C = unique(B,'rows','stable');
#fid = fopen('Stest.txt','w');
#fprintf(fid,'%s\n',C');
#fclose(fid);
#EOF

while read h V2
do
cd "${h:0:33}"
echo ${h:0:33}
for i in 2*PE01*H*
do
if [[ "${#i}" -gt 5 ]]; then
echo $i
wav_year=${i:0:4}
wav_jday=${i:4:3}
jtemp=`echo "${wav_jday} - 1" | bc -l`
wav_md=`date -d "${wav_year}-01-01 ${jtemp} days" +%m/%d`
wav_hwr=${i:7:2}
wav_min=${i:9:2}
wav_sec=${i:11:2}
wav_msec=${i:14:2}
echo ${wav_year}/${wav_md}T${wav_hwr}:${wav_min}:${wav_sec}.${wav_msec}
temp=${i:(-8):4}
stnm=${temp//.}
comp=${i:(-3):3}
pzf=`find /home/user/env_EGEL/EGEL_PZs/* -name "SAC*${stnm}_${comp}*"` 
j=0
for j in ${pzf};
do
A_j=${j: -18:3}
A_y=${j: -23:4}
A_h=${j: -14:2}
A_m=${j: -11:2}
A_s=${j: -8:2}
A_j=$(echo "$A_j-1" | bc)
A_s=$(echo "$A_s-1" | bc)
A_m=$(echo "$A_m-1" | bc)
A_h=$(echo "$A_h-1" | bc)
A_us=${j: -5:5}
A_md=`date -d "${A_y}-01-01 ${A_j} days" +%m/%d`
A_pzf=`date -d "${A_y}/${A_md}T${A_h}:${A_m}:${A_s}" +%s`

B_j=${j: -41:3}
B_y=${j: -46:4}
B_h="00"
B_m="00"
B_s="00"
B_us="00000"
B_j=$(echo "$B_j-1" | bc)
B_md=`date -d "${B_y}-01-01 ${B_j} days" +%m/%d`
B_pzf=`date -d "${B_y}/${B_md}T${B_h}:${B_m}:${B_s}.${B_us}" +%s`

A_o=`date -d "${wav_year}/${wav_md}T${wav_hwr}:${wav_min}:${wav_sec}.${wav_msec}" +%s`
#echo $A_o $A_pzf $B_pzf
if [[ "${A_o}" -ge "${B_pzf}" ]] && [[ "${A_o}" -le "${A_pzf}" ]]; then
o_var=$j
fi
done
pzf=$o_var
echo $pzf
if [[ "${comp:0:1}" = "B" ]]; then
	c1="0.1"
	c2="1"
	c3="8.5"
	c4="9.5"
elif [[ "${comp:0:1}" = "S" ]]; then
	c1="0.1"
	c2="1"
	c3="18"
	c4="22"
else
	c1="0.1"
	c2="1"
	c3="40"
	c4="48"
fi

rm -f ${stnm}_${comp}_1-2.sacii ${stnm}_${comp}_1-2.ascii
rm -f ${stnm}_${comp}_2-4.sacii ${stnm}_${comp}_2-4.ascii
rm -f ${stnm}_${comp}_4-8.sacii ${stnm}_${comp}_4-8.ascii
rm -f ${stnm}_${comp}_8-16.sacii ${stnm}_${comp}_8-16.ascii
rm -f ${stnm}_${comp}_16-32.sacii ${stnm}_${comp}_16-32.ascii
rm -f ${stnm}_${comp}.sacii ${stnm}_${comp}.ascii

sac <<EOF
   wild echo off
   read $i
   ch lovrok true
   r
   rmean
   rtrend
   trans from polezero s ${pzf} to vel freqlimits ${c1} ${c2} ${c3} ${c4}
   write SAC N_${stnm}_${comp}.sacii
   quit
EOF

if [[ "${comp:0:1}" = "B" ]]; then
 sac <<KOD
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 1 2
w SAC N_${stnm}_${comp}_1-2.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 2 4
w SAC N_${stnm}_${comp}_2-4.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 4 8
w SAC N_${stnm}_${comp}_4-8.sacii
quit
KOD
elif [[ "${comp:0:1}" = "S" ]]; then
sac <<MOD
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 1 2
w SAC N_${stnm}_${comp}_1-2.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 2 4
w SAC N_${stnm}_${comp}_2-4.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 4 8
w SAC N_${stnm}_${comp}_4-8.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 8 16
w SAC N_${stnm}_${comp}_8-16.sacii
quit
MOD
else
 sac <<SOD
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 1 2
w SAC N_${stnm}_${comp}_1-2.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 2 4
w SAC N_${stnm}_${comp}_2-4.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 4 8
w SAC N_${stnm}_${comp}_4-8.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 8 16
w SAC N_${stnm}_${comp}_8-16.sacii
r N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 16 32
w SAC N_${stnm}_${comp}_16-32.sacii
quit
SOD
fi
fi
done
cd ../..
done < Stest.txt
