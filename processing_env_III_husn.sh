#!/bin/bash

#rm -f ~/env_EGEL/intgrl_husn_compil_fin2.txt

for h in /home/user/env_EGEL/husn_events/2019-02-{05..28..1}* /home/user/env_EGEL/husn_events/2019-{03..05..1}*
do
cd "$h"
echo $h
#rm -f *.ascii
rm -f outful*.txt
#rm -f env_integral.txt
cw=`grep -A1 "Date" ./husn*.pick | tail -n1`
Oymd=`echo ${cw} | awk '{print $1}'`
Ohms=`echo ${cw} | awk '{print $2}'`
Odep=`echo ${cw} | awk '{print $9}'`
Oyear=${Oymd:0:4}
Omnth=${Oymd:5:2}
Odate=${Oymd:8:2}
Ojday=`date -d "${Oyear}/${Omnth}/${Odate}" +%j`
Ohr=${Ohms:0:2}
Omin=${Ohms:3:2}
Osec=${Ohms:6:2}
Omsec=${Osms:9:2}
Omsec=${Omsec}"0"
for i in N_*_*HZ_*.sacii 
do
echo $i
if [[ "${#i}" -gt 16 ]]; then
file2=${i/HZ_/HN_}
file3=${i/HZ_/HE_}
#echo $i $file2 $file3
temp1=${i:(-11):12}
temp2=${temp1//.sacii}
temp3=${temp2//Z}
frange=${temp3//_}
#echo $frange
stnm=${i:2:4}
stnm=${stnm//_}

Shrminss=`grep " ${stnm} " ./husn*.pick | grep " S_" | awk '{ print $5 }'`
Styp=`grep " ${stnm} " ./husn*.pick | grep " S_" | awk '{ print $4 }'`
#echo $Shrminss
if [[ "${#Shrminss}" -gt "0" ]]; then
Syear=${Oyear}
Smnth=${Omnth}
Sdate=${Odate}
Sjday=`date -d "${Syear}/${Smnth}/${Sdate}" +%j`
Shr=${Shrminss:0:2}
Smin=${Shrminss:3:2}
Ssec=${Shrminss:6:2}
Smsec=${Shrminss:9:3}

absStime=`date -d "${Syear}/${Smnth}/${Sdate}T${Shr}:${Smin}:${Ssec}.${Smsec}" +%s.%N`

Phrminss=`grep " ${stnm} " ./husn*.pick | grep " P_" | awk '{ print $5 }'`
Pyear=${Oyear}
Pmnth=${Omnth}
Pdate=${Odate}
Pjday=`date -d "${Pyear}/${Pmnth}/${Pdate}" +%j`
Phr=${Phrminss:0:2}
Pmin=${Phrminss:3:2}
Psec=${Phrminss:6:2}
Pmsec=${Phrminss:9:3}

if [[ "${#Phrminss}" -eq "0" ]]; then
Phr=${Ohms:0:2}
Pmin=${Ohms:3:2}
Psec=${Ohms:6:2}
Pmsec=${Ohms:9:2}
Pmsec=${Pmsec}"0"
fi

absPtime=`date -d "${Pyear}/${Pmnth}/${Pdate}T${Phr}:${Pmin}:${Psec}.${Pmsec}" +%s.%N`

#echo $Oyear $Omnth $Odate $Ojday $Ohr $Omin $Osec
#echo $Syear $Smnth $Sdate $Sjday $Shr $Smin $Ssec $Smsec

Sdist=`grep " ${stnm} " ./husn*.pick | grep " S_" | awk '{ print $2 }'`
Sdist=$(echo "${Sdist}*111.1949" | bc -l)
Sdist=$(echo "scale=4;sqrt(${Sdist}^2 + ${Odep}^2)" | bc -l)
#echo $Sdist
echo $stnm $Sdist >> st_dist.txt

sac <<EOF
r ${i}
ch lovrok true
wh
ch A GMT ${Syear} ${Sjday} ${Shr} ${Smin} ${Ssec} ${Smsec}
wh
ch O GMT ${Oyear} ${Ojday} ${Ohr} ${Omin} ${Osec} ${Omsec}
wh
r ${file2}
ch lovrok true
wh
ch A GMT ${Syear} ${Sjday} ${Shr} ${Smin} ${Ssec} ${Smsec}
wh
ch O GMT ${Oyear} ${Ojday} ${Ohr} ${Omin} ${Osec} ${Omsec}
wh
r ${file3}
ch lovrok true
wh
ch A GMT ${Syear} ${Sjday} ${Shr} ${Smin} ${Ssec} ${Smsec}
wh
ch O GMT ${Oyear} ${Ojday} ${Ohr} ${Omin} ${Osec} ${Omsec}
wh
quit
EOF

diff=$(echo "$absStime - $absPtime + 20" | bc -l)
#echo $diff

FSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${i} | awk '{ print $2 }'`
FSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file2}| awk '{ print $2 }'`
FSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file3}| awk '{ print $2 }'`
#FSNR_Z2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
#FSNR_N2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
#FSNR_E2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
#FSNR_Z3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
#FSNR_N3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
#FSNR_E3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`

diff=$(echo "$absStime - $absPtime + 40" | bc -l)
#echo $diff

SSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${i} | awk '{ print $2 }'`
SSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file2}| awk '{ print $2 }'`
SSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file3}| awk '{ print $2 }'`
#SSNR_Z2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
#SSNR_N2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
#SSNR_E2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
#SSNR_Z3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
#SSNR_N3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
#SSNR_E3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`

diff=$(echo "$absStime - $absPtime + 60" | bc -l)
#echo $diff

TSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${i} | awk '{ print $2 }'`
TSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file2}| awk '{ print $2 }'`
TSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file3}| awk '{ print $2 }'`
#TSNR_Z2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
#TSNR_N2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
#TSNR_E2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
#TSNR_Z3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
#TSNR_N3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
#TSNR_E3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`


SNR_Z1=`perl -e "use List::Util qw[min max];print max(${FSNR_Z1},${SSNR_Z1},${TSNR_Z1})"`
SNR_N1=`perl -e "use List::Util qw[min max];print max(${FSNR_N1},${SSNR_N1},${TSNR_N1})"`
SNR_E1=`perl -e "use List::Util qw[min max];print max(${FSNR_E1},${SSNR_E1},${TSNR_E1})"`
#SNR_Z2=`perl -e "use List::Util qw[min max];print max(${FSNR_Z2},${SSNR_Z2},${TSNR_Z2})"`
#SNR_N2=`perl -e "use List::Util qw[min max];print max(${FSNR_N2},${SSNR_N2},${TSNR_N2})"`
#SNR_E2=`perl -e "use List::Util qw[min max];print max(${FSNR_E2},${SSNR_E2},${TSNR_E2})"`
#SNR_Z3=`perl -e "use List::Util qw[min max];print max(${FSNR_Z3},${SSNR_Z3},${TSNR_Z3})"`
#SNR_N3=`perl -e "use List::Util qw[min max];print max(${FSNR_N3},${SSNR_N3},${TSNR_N3})"`
#SNR_E3=`perl -e "use List::Util qw[min max];print max(${FSNR_E3},${SSNR_E3},${TSNR_E3})"`

if (( $(echo "$SNR_Z1 > 3" | bc -l) )) && (( $(echo "$SNR_N1 > 3" | bc -l) )) && (( $(echo "$SNR_E1 > 3" | bc -l) )); then
#&& (( $(echo "$SNR_Z2 > 3" | bc -l) )) && (( $(echo "$SNR_Z3 > 3" | bc -l) )) 

#&& (( $(echo "$SNR_N2 > 3" | bc -l) )) && (( $(echo "$SNR_N3 > 3" | bc -l) )) 

#&& (( $(echo "$SNR_E2 > 3" | bc -l) )) && (( $(echo "$SNR_E3 > 3" | bc -l) )); then

  echo success

  if [[ "${frange}" = "1-2" ]]; then
    delf="1"
  elif [[ "${frange}" = "2-4" ]]; then
    delf="2"
  elif [[ "${frange}" = "4-8" ]]; then
    delf="4"
  elif [[ "${frange}" = "8-16" ]]; then
    delf="8"
  else
    delf="16"
  fi
  rm -f ${i/HZ/env}
sac <<EOF
r ${i}
hilbert
sqr
write sac ${i/N_/hilsqr_}
r ${i}
sqr
write sac ${i/N_/sqr_}
r ${file2}
hilbert
sqr
write sac ${file2/N_/hilsqr_}
r ${file2}
sqr
write sac ${file2/N_/sqr_}
r ${file3}
hilbert
sqr
write sac ${file3/N_/hilsqr_}
r ${file3}
sqr
write sac ${file3/N_/sqr_}
r ${i/N_/hilsqr_}
addf ${i/N_/sqr_} ${file2/N_/hilsqr_} ${file2/N_/sqr_} ${file3/N_/hilsqr_} ${file3/N_/sqr_}
div 8
mul 2.83
div ${delf}
write sac ${i/HZ/env}
r ${i/HZ/env}
ch lovrok true
wh
ch T0 GMT ${Pyear} ${Pjday} ${Phr} ${Pmin} ${Psec} ${Pmsec}
wh
cut T0 -50 -45
r
write alpha cut1.ascii
cut off
r ${i/HZ/env}
cut T0 -45 -40
r
write alpha cut2.ascii
cut off
r ${i/HZ/env}
cut T0 -40 -35
r
write alpha cut3.ascii
cut off
r ${i/HZ/env}
cut T0 -35 -30
r
write alpha cut4.ascii
cut off
r ${i/HZ/env}
cut T0 -30 -25
r
write alpha cut5.ascii
cut off
r ${i/HZ/env}
cut T0 -25 -20
r
write alpha cut6.ascii
cut off
r ${i/HZ/env}
cut T0 -20 -15
r
write alpha cut7.ascii
cut off
r ${i/HZ/env}
cut T0 -15 -10
r
write alpha cut8.ascii
cut off
quit
EOF
rm -f ${i/N_/hilsqr_} ${i/N_/sqr_} ${file2/N_/hilsqr_} ${file2/N_/sqr_} ${file3/N_/hilsqr_} ${file3/N_/sqr_}
# ${i} ${file2} ${file3}
avg1=`head -n12 cut1.ascii | tail -n1 | awk '{ print $2}'`
avg2=`head -n12 cut2.ascii | tail -n1 | awk '{ print $2}'`
avg3=`head -n12 cut3.ascii | tail -n1 | awk '{ print $2}'`
avg4=`head -n12 cut4.ascii | tail -n1 | awk '{ print $2}'`
avg5=`head -n12 cut5.ascii | tail -n1 | awk '{ print $2}'`
avg6=`head -n12 cut6.ascii | tail -n1 | awk '{ print $2}'`
avg7=`head -n12 cut7.ascii | tail -n1 | awk '{ print $2}'`
avg8=`head -n12 cut8.ascii | tail -n1 | awk '{ print $2}'`
avg_min=`perl -e "use List::Util qw[min max];print min(${avg1},${avg2},${avg3},${avg4},${avg5},${avg6},${avg7},${avg8})"`
#echo $avg1 $avg2 $avg3 $avg4 $avg_min
avg1=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg1"`
avg2=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg2"`
avg3=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg3"`
avg4=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg4"`
avg5=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg5"`
avg6=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg6"`
avg7=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg7"`
avg8=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg8"`
avg_min=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg_min"`

if (( $(echo "${avg1} == ${avg_min}" | bc -l) )); then
 #echo ${avg1}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 50" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
elif (( $(echo "${avg2} == ${avg_min}" | bc -l) )); then
 #echo ${avg2}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 45" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
elif (( $(echo "${avg3} == ${avg_min}" | bc -l) )); then
 #echo ${avg3}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 40" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
elif (( $(echo "${avg4} == ${avg_min}" | bc -l) )); then
 #echo ${avg4}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 35" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
elif (( $(echo "${avg5} == ${avg_min}" | bc -l) )); then
 #echo ${avg5}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 30" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
elif (( $(echo "${avg6} == ${avg_min}" | bc -l) )); then
 #echo ${avg6}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 25" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
elif (( $(echo "${avg7} == ${avg_min}" | bc -l) )); then
# echo ${avg7}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 20" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
else
# echo ${avg8}
 for (( j=1; j<=19; j++ ))
 do
   diff=$(echo "$absStime - $absPtime + 15" | bc -l)
   Sstrt=`echo "${j-1}*5+10" | bc`
   SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
   if (( $(echo "$SNR_env < 3" | bc -l) )); then
     cut_W=${Sstrt}
     break;
   else
     cut_W=${Sstrt}
   fi
 done
fi
rm -f cut*.ascii
rm -f ${i/HZ/Swin} ${i/HZ/Codawin} SNR_coda.sacii ${i/HZ/env2} SNR_coda.ascii out_SNR_coda.txt
temp=${i/HZ/envwin}
sac <<EOF
r ${i/HZ/env}
ch T0 GMT ${Syear} ${Sjday} ${Shr} ${Smin} ${Ssec} ${Smsec}
wh
cut T0 -1 ${cut_W}
r
write sac ${i/HZ/envful}
write alpha ${temp/sacii/ascii}
cut off
quit
EOF
tail -n +31 ${temp/sacii/ascii} > test.txt
awk '{for (i=1;i<=NF;i++) print $i}' test.txt >> outful_${temp/sacii/txt}
echo $Sdist >> outful_${temp/sacii/txt}

#: <<SOD

#SOD
rm -f ${temp/sacii/ascii} ${i/HZ/env}
fi
fi
fi
#rm -f ${i} ${file2} ${file3}

done
if [[ "${#Shrminss}" -gt "0" ]]; then
octave <<EOF
%profile on;
alpha = dir("outful_N*");
beta = struct2cell(alpha);
for p = 1:length(alpha)
A = load(beta{(p-1)*6+1});
b = bartlett(101);
D = conv(A(1:end-1),b/sum(b),'same');
pkg load signal
[q,loc] = findpeaks(D(1001:end));
E = 1./D;
[qd,locd] = findpeaks(E(1001:end));
loc = loc + 1000;
locd = locd + 1000;
%length(loc)
%length(locd)
x = 1:1:length(D);
%%plot(x,D); hold on; plot(x(loc),q,'ro*'); plot(x(locd),D(locd),'o'); hold off;
jmp = zeros(length(loc),1);
idx = zeros(length(loc),1);
for i=1:length(loc)
	sw = find(locd<loc(i));
	jmp(i) = 0;
        idx(i) = 0;
	if (length(sw) > 0)
		idx(i) = locd(sw(end));
		jmp(i) = D(loc(i))/D(idx(i));
	end
end
rw = find(jmp > 4);
cut_p = length(D) - 1;
if (length(rw) > 0)
	cut_p = idx(rw(1));
end

fid = fopen(beta{(p-1)*6+1},'w');
fprintf(fid,"%e\n",D(1:cut_p)');
fprintf(fid,"%f",A(end));
fclose(fid);
end
%profile off;
%T = profile("info");
%profshow(T);
EOF
fi
rm -f N_*HH*.sacii
#touch env_integral.txt
#SOD
#cat env_integral.txt >> /home/user/env_EGEL/intgrl_husn_compil_fin2.txt
cd ..
done
