#!/bin/bash

rm -f ~/env_EGEL/intgrl_egel_compil_sc_nr.txt

for h in /home/user/env_EGEL/EGEL_events/*
do
cd $h
echo $h
rm -f *.ascii
rm -f out*.txt
rm -f env_integral.txt

Oyear=`grep "GEOGRAPHIC" *.hyp | awk '{ print $3 }'`
Omnth=`grep "GEOGRAPHIC" *.hyp | awk '{ print $4 }'`
Odate=`grep "GEOGRAPHIC" *.hyp | awk '{ print $5 }'`
Ojday=`date -d "${Oyear}/${Omnth}/${Odate}" +%j`
Ohr=`grep "GEOGRAPHIC" *.hyp | awk '{ print $6 }'`
Omin=`grep "GEOGRAPHIC" *.hyp | awk '{ print $7 }'`
Osms=`grep "GEOGRAPHIC" *.hyp | awk '{ print $8 }'`
Odep=`grep "GEOGRAPHIC" *.hyp | awk '{ print $14 }'`
Osec=${Osms:0:2}
Osec=${Osec//.}
Omsec=$(echo "${Osms//.}%1000" | bc)

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

Symd=`grep "^${stnm}" ./egel.2*.hyp | grep "? S" | awk '{ print $7 }'`
Shrmin=`grep "^${stnm}" ./egel.2*.hyp | grep "? S" | awk '{ print $8 }'`
Ssms=`grep "^${stnm}" ./egel.2*.hyp | grep "? S" | awk '{ print $9 }'`
Syear=${Symd:0:4}
Smnth=${Symd:4:2}
Sdate=${Symd:6:2}
Sjday=`date -d "${Syear}/${Smnth}/${Sdate}" +%j`
Shr=${Shrmin:0:2}
Smin=${Shrmin:2:2}
Ssec=${Ssms:0:2}
Ssec=${Ssec//.}
Smsec=${Ssms:(-4):3}

absStime=`date -d "${Syear}/${Smnth}/${Sdate}T${Shr}:${Smin}:${Ssec}.${Smsec}" +%s.%N`

Pymd=`grep "^${stnm}" ./egel.2*.hyp | grep "? P" | awk '{ print $7 }'`
Phrmin=`grep "^${stnm}" ./egel.2*.hyp | grep "? P" | awk '{ print $8 }'`
Psms=`grep "^${stnm}" ./egel.2*.hyp | grep "? P" | awk '{ print $9 }'`
Pyear=${Pymd:0:4}
Pmnth=${Pymd:4:2}
Pdate=${Pymd:6:2}
Pjday=`date -d "${Pyear}/${Pmnth}/${Pdate}" +%j`
Phr=${Phrmin:0:2}
Pmin=${Phrmin:2:2}
Psec=${Psms:0:2}
Psec=${Psec//.}
Pmsec=${Psms:(-4):3}

if [[ "${#Pymd}" -eq "0" ]]; then
Pymd=${Symd}
Phrmin=$Shrmin
Psms=$Ssms
Pyear=${Pymd:0:4}
Pmnth=${Pymd:4:2}
Pdate=${Pymd:6:2}
Pjday=`date -d "${Pyear}/${Pmnth}/${Pdate}" +%j`
Phr=${Phrmin:0:2}
Pmin=${Phrmin:2:2}
Psec=${Psms:0:2}
Psec=${Psec//.}
Pmsec=${Psms:(-4):3}
fi

absPtime=`date -d "${Pyear}/${Pmnth}/${Pdate}T${Phr}:${Pmin}:${Psec}.${Pmsec}" +%s.%N`

#echo $Oyear $Omnth $Odate $Ojday $Ohr $Omin $Osec
#echo $Syear $Smnth $Sdate $Sjday $Shr $Smin $Ssec $Smsec

Sdist=`grep "^${stnm}" ./egel.2*.hyp | grep "? S" | awk '{ print $22 }'`
Sdist=$(echo "scale=4;sqrt(${Sdist}^2 + ${Odep}^2)" | bc -l)
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

FSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${i} | awk '{ print $2 }'`
FSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
FSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
FSNR_Z2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
FSNR_N2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
FSNR_E2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
FSNR_Z3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
FSNR_N3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
FSNR_E3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
FSNR_Z4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
FSNR_N4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
FSNR_E4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`

diff=$(echo "$absStime - $absPtime + 40" | bc -l)
#echo $diff

SSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${i} | awk '{ print $2 }'`
SSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
SSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
SSNR_Z2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
SSNR_N2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
SSNR_E2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
SSNR_Z3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
SSNR_N3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
SSNR_E3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
SSNR_Z4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
SSNR_N4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
SSNR_E4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`

diff=$(echo "$absStime - $absPtime + 60" | bc -l)
#echo $diff

TSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${i} | awk '{ print $2 }'`
TSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
TSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
TSNR_Z2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
TSNR_N2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
TSNR_E2=`sac_snr -A-2 -S15 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
TSNR_Z3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
TSNR_N3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
TSNR_E3=`sac_snr -A-2 -S30 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`
TSNR_Z4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${i}| awk '{ print $2 }'`
TSNR_N4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${file2}| awk '{ print $2 }'`
TSNR_E4=`sac_snr -A-2 -S45 -N-${diff} -W15 ${file3}| awk '{ print $2 }'`


SNR_Z1=`perl -e "use List::Util qw[min max];print max(${FSNR_Z1},${SSNR_Z1},${TSNR_Z1})"`
SNR_N1=`perl -e "use List::Util qw[min max];print max(${FSNR_N1},${SSNR_N1},${TSNR_N1})"`
SNR_E1=`perl -e "use List::Util qw[min max];print max(${FSNR_E1},${SSNR_E1},${TSNR_E1})"`
SNR_Z2=`perl -e "use List::Util qw[min max];print max(${FSNR_Z2},${SSNR_Z2},${TSNR_Z2})"`
SNR_N2=`perl -e "use List::Util qw[min max];print max(${FSNR_N2},${SSNR_N2},${TSNR_N2})"`
SNR_E2=`perl -e "use List::Util qw[min max];print max(${FSNR_E2},${SSNR_E2},${TSNR_E2})"`
SNR_Z3=`perl -e "use List::Util qw[min max];print max(${FSNR_Z3},${SSNR_Z3},${TSNR_Z3})"`
SNR_N3=`perl -e "use List::Util qw[min max];print max(${FSNR_N3},${SSNR_N3},${TSNR_N3})"`
SNR_E3=`perl -e "use List::Util qw[min max];print max(${FSNR_E3},${SSNR_E3},${TSNR_E3})"`
SNR_Z4=`perl -e "use List::Util qw[min max];print max(${FSNR_Z4},${SSNR_Z4},${TSNR_Z4})"`
SNR_N4=`perl -e "use List::Util qw[min max];print max(${FSNR_N4},${SSNR_N4},${TSNR_N4})"`
SNR_E4=`perl -e "use List::Util qw[min max];print max(${FSNR_E4},${SSNR_E4},${TSNR_E4})"`

if (( $(echo "$SNR_Z1 > 3" | bc -l) )) && (( $(echo "$SNR_Z2 > 2" | bc -l) )) && (( $(echo "$SNR_Z3 > 2" | bc -l) )) && (( $(echo "$SNR_Z4 > 1" | bc -l) )) && (( $(echo "$SNR_N1 > 3" | bc -l) )) && (( $(echo "$SNR_N2 > 2" | bc -l) )) && (( $(echo "$SNR_N3 > 2" | bc -l) )) && (( $(echo "$SNR_N4 > 1" | bc -l) )) && (( $(echo "$SNR_E1 > 3" | bc -l) )) && (( $(echo "$SNR_E2 > 2" | bc -l) )) && (( $(echo "$SNR_E3 > 2" | bc -l) )) && (( $(echo "$SNR_E4 > 1" | bc -l) )); then

SNR_F1=`echo "scale=4;(${SNR_Z1} + ${SNR_E1} + ${SNR_N1})/3" | bc -l`
SNR_F2=`echo "scale=4;(${SNR_Z2} + ${SNR_E2} + ${SNR_N2})/3" | bc -l`
SNR_F3=`echo "scale=4;(${SNR_Z3} + ${SNR_E3} + ${SNR_N3})/3" | bc -l`

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
quit
EOF
rm -f ${i/N_/hilsqr_} ${i/N_/sqr_} ${file2/N_/hilsqr_} ${file2/N_/sqr_} ${file3/N_/hilsqr_} ${file3/N_/sqr_}
#  ${i} ${file2} ${file3}
: <<SOD

sac <<EOF
r ${i/HZ/env}
ch lovrok true
wh
ch T0 GMT ${Pyear} ${Pjday} ${Phr} ${Pmin} ${Psec} ${Pmsec}
wh
cut T0 -50 -40
r
write alpha cut1.ascii
cut off
r ${i/HZ/env}
cut T0 -40 -30
r
write alpha cut2.ascii
cut off
r ${i/HZ/env}
cut T0 -30 -20
r
write alpha cut3.ascii
cut off
r ${i/HZ/env}
cut T0 -20 -10
r
write alpha cut4.ascii
cut off
quit
EOF
#avg1=`head -n12 cut1.ascii | tail -n1 | awk '{ print $2}'`
#avg2=`head -n12 cut2.ascii | tail -n1 | awk '{ print $2}'`
#avg3=`head -n12 cut3.ascii | tail -n1 | awk '{ print $2}'`
#avg4=`head -n12 cut4.ascii | tail -n1 | awk '{ print $2}'`
#avg_min=`perl -e "use List::Util qw[min max];print min(${avg1},${avg2},${avg3},${avg4})"`
#echo $avg1 $avg2 $avg3 $avg4 $avg_min
#rm -f cut*.ascii
#sac <<EOF
#r ${i/HZ/env}
#sub ${avg_min}
#write sac ${i/HZ/env}
#quit
#EOF
SOD
rm -f ${i/HZ/Swin} ${i/HZ/Codawin} SNR_coda.sacii ${i/HZ/env2} SNR_coda.ascii out_SNR_coda.txt
temp=${i/HZ/envwin}
sac <<EOF
r ${i/HZ/env}
ch T0 GMT ${Syear} ${Sjday} ${Shr} ${Smin} ${Ssec} ${Smsec}
wh
cut T0 -1 100
r
write sac ${i/HZ/envwin}
write alpha ${temp/sacii/ascii}
cut off
quit
EOF
tail -n +31 ${temp/sacii/ascii} > test.txt
awk '{for (i=1;i<=NF;i++) print $i}' test.txt >> out_${temp/sacii/txt}
echo $Sdist >> out_${temp/sacii/txt}
echo $SNR_F1 >> out_${temp/sacii/txt}
echo $SNR_F2 >> out_${temp/sacii/txt}
echo $SNR_F3 >> out_${temp/sacii/txt}
rm -f ${temp/sacii/ascii} ${i/HZ/env}
fi
fi
done


octave <<EOF

clc;
close all;

list = dir('out_N*envwin*.txt');

filnam = sprintf('env_integral.txt');
fid = fopen(filnam,'w');

for i_iter = 1:length(list)
    
ay = strrep(list(i_iter).name,'out_N_','');

stnm = ay(1:4);
bnd = ay(end-6:end-4);
if (ay(4)=='_')
    stnm = ay(1:3);
end
if (bnd(1)=='-')
   bnd = ay(end-7:end-4);
   if (bnd(1)=='6')
	bnd = ay(end-8:end-4);
   end
end

temp = load(list(i_iter).name);

sdist = num2str(temp(end-3))

snr1 = num2str(temp(end-2));

snr2 = num2str(temp(end-1));

snr3 = num2str(temp(end));

dta = temp(1:end-4);

if length(dta) > 10000
  fsamp = 100;
elseif (length(dta) > 5000) && (length(dta) < 10000)
  fsamp = 50;
  x = linspace(-1,100,length(dta));
  x_nw = linspace(-1,100,10100);
  dta_nw = interp1(x,dta,x_nw);
  dta = dta_nw;
else
  fsamp = 20;
  x = linspace(-1,100,length(dta));
  x_nw = linspace(-1,100,10100);
  dta_nw = interp1(x,dta,x_nw);
  dta = dta_nw;
end

fsamp=100;

U4 = trapz(dta(fsamp*51:fsamp*61))/(10*fsamp);
U1 = trapz(dta(fsamp*1:fsamp*16))/U4;
U2 = trapz(dta(fsamp*16+1:fsamp*31))/U4;
U3 = trapz(dta(fsamp*31+1:fsamp*46))/U4;

U1 = num2str(U1);
U2 = num2str(U2);
U3 = num2str(U3);

AB = {stnm, sdist, bnd, 'U1_norm', U1, snr1; stnm, sdist, bnd, 'U2_norm', U2, snr2; stnm, sdist, bnd, 'U3_norm', U3, snr3};

fprintf(fid,"%s %s %s %s %s %s\n",AB'{:});

end
fclose(fid);
close all;
clear all;
quit
EOF
#SOD
cat env_integral.txt >> /home/user/env_EGEL/intgrl_egel_compil_sc_nr.txt
cd ..
done
