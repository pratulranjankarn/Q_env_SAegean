#!/bin/bash

OLDIFS=${IFS}

#: <<SOD

for h in /home/user/env_EGEL/kef_events/20* #/home/user/env_EGEL/cornoth_events/20* 
#/home/user/env_EGEL/husn_events/2012-{01..12..1}* /home/user/env_EGEL/husn_events/{2013..2018..1}* /home/user/env_EGEL/husn_events/2019-{01..05..1}*
##/home/user/env_EGEL/husn_events/2012-10-06T04M00S1* 
##/home/user/env_EGEL/husn_events/2012-{10..12..1}* /home/user/env_EGEL/husn_events/{2013..2018..1}* /home/user/env_EGEL/husn_events/2019-{01..05..1}*
##/home/user/env_EGEL/husn_events/2013-01-01T01M18S3* 
do
##	cd $h
	echo $h	
	rm -f ${h}/outful*.txt
	rm -f ${h}/*.ascii
	hyp_f=`find ${h} -name "husn*.hyp"`
##	echo $hyp_f
	if [ "${#hyp_f}" -gt 0 ]; then

		O_info=`grep "GEOGRAPHIC" ${hyp_f}`
		read -ra O_arr <<< "${O_info}"
##		echo ${O_arr[0]} ${O_arr[1]} ${O_arr[2]} ${O_arr[3]} ${O_arr[4]} ${O_arr[5]} ${O_arr[6]} ${O_arr[7]}
		Oyear=${O_arr[2]}
		Omnth=${O_arr[3]}
		Odate=${O_arr[4]}
##		echo D0
		Ojday=`date -d "${Oyear}/${Omnth}/${Odate}" +%j`
		Ohr=${O_arr[5]}
		Omin=${O_arr[6]}
		Osms=${O_arr[7]}
		Odep=${O_arr[13]}
		Osec=${Osms:0:2}
		Osec=${Osec//.}
		Omsec=$(echo "${Osms//.}%1000" | bc)


		for i in ${h}/N_*_*HZ_*.sacii 
		do
			file2=${i/HZ_/HN_}
			file3=${i/HZ_/HE_}
			if [ "${#i}" -gt 70 ] && [ -a ${file2} ] && [ -a ${file3} ]; then
				echo $i 
##				$file2 $file3
				IFS="_"
				read -ra f_info <<< "${i}"
				IFS=${OLDIFS}
				temp1=${f_info[5]}
				frange=${temp1//.sacii}
##				echo $frange
				stnm=${f_info[3]}
				comp=${f_info[4]}

				S_info=`grep "^${stnm}" ${hyp_f} | grep "? S"`
				read -ra S_arr <<< "${S_info}"
##				echo ${S_arr[0]} ${S_arr[1]} ${S_arr[2]} ${S_arr[3]} ${S_arr[4]} ${S_arr[5]} ${S_arr[6]} ${S_arr[7]} ${S_arr[8]} ${S_arr[9]} ${S_arr[10]} ${S_arr[11]}

				if [ "${#S_arr}" -gt 0 ]; then
##					echo D1
					Syear=${S_arr[6]:0:4}
					Smnth=${S_arr[6]:4:2}
					Sdate=${S_arr[6]:6:2}
					Sjday=`date -d "${Syear}/${Smnth}/${Sdate}" +%j`
					Shr=${S_arr[7]:0:2}
					Smin=${S_arr[7]:2:2}
					Ssec=${S_arr[8]:0:2}
					Ssec=${Ssec//.}
					Smsec=${S_arr[8]:(-4):3}
					
##					echo D2
					absStime=`date -d "${Syear}/${Smnth}/${Sdate}T${Shr}:${Smin}:${Ssec}.${Smsec}" +%s.%N`
					
					P_info=`grep "^${stnm}" ${hyp_f} | grep "? P"`
					if [[ "${#P_info}" -gt "0" ]]; then
						read -ra P_arr <<< "${P_info}"
						Pyear=${P_arr[6]:0:4}
						Pmnth=${P_arr[6]:4:2}
						Pdate=${P_arr[6]:6:2}
##						echo D3
						Pjday=`date -d "${Pyear}/${Pmnth}/${Pdate}" +%j`
						Phr=${P_arr[7]:0:2}
						Pmin=${P_arr[7]:2:2}
						Psec=${P_arr[8]:0:2}
						Psec=${Psec//.}
						Pmsec=${P_arr[8]:(-4):3}
					else
						Pyear=${Syear}
						Pmnth=${Smnth}
						Pdate=${Sdate}
##						echo D4
						Pjday=`date -d "${Pyear}/${Pmnth}/${Pdate}" +%j`
						Phr=${Shr}
						Pmin=${Smin}
						Psec=${Ssec}
						Pmsec=${Smsec}
					fi
					
##					echo D5
					absPtime=`date -d "${Pyear}/${Pmnth}/${Pdate}T${Phr}:${Pmin}:${Psec}.${Pmsec}" +%s.%N`
					
##					echo $Oyear $Omnth $Odate $Ojday $Ohr $Omin $Osec $Omsec
##					echo $Syear $Smnth $Sdate $Sjday $Shr $Smin $Ssec $Smsec

					Sdist=${S_arr[21]}
					if [ "${Sdist}" = "0.0000" ]; then
						Sdist="400"
					fi
					Sdist=$(echo "scale=4;sqrt(${Sdist}^2 + ${Odep}^2)" | bc -l)
					rm -f st_dist.txt


sac <<EOF
rh ${i} ${file2} ${file3}
ch lovrok true
ch A GMT ${Syear} ${Sjday} ${Shr} ${Smin} ${Ssec} ${Smsec}
ch O GMT ${Oyear} ${Ojday} ${Ohr} ${Omin} ${Osec} ${Omsec}
wh
quit
EOF

					diff=$(echo "$absStime - $absPtime + 20" | bc -l)
##					echo $diff

					FSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${i} | awk '{ print $2 }'`
					FSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file2}| awk '{ print $2 }'`
					FSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file3}| awk '{ print $2 }'`
					
					diff=$(echo "$absStime - $absPtime + 40" | bc -l)
##					echo $diff
					
					SSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${i} | awk '{ print $2 }'`
					SSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file2}| awk '{ print $2 }'`
					SSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file3}| awk '{ print $2 }'`
					
					diff=$(echo "$absStime - $absPtime + 60" | bc -l)
##					echo $diff
					
					TSNR_Z1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${i} | awk '{ print $2 }'`
					TSNR_N1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file2}| awk '{ print $2 }'`
					TSNR_E1=`sac_snr -A-2 -S-1 -N-${diff} -W10 ${file3}| awk '{ print $2 }'`
					
					
					SNR_Z1=`perl -e "use List::Util qw[min max];print max(${FSNR_Z1},${SSNR_Z1},${TSNR_Z1})"`
					SNR_N1=`perl -e "use List::Util qw[min max];print max(${FSNR_N1},${SSNR_N1},${TSNR_N1})"`
					SNR_E1=`perl -e "use List::Util qw[min max];print max(${FSNR_E1},${SSNR_E1},${TSNR_E1})"`
					
					if (( $(echo "$SNR_Z1 > 3" | bc -l) )) && (( $(echo "$SNR_N1 > 3" | bc -l) )) && (( $(echo "$SNR_E1 > 3" | bc -l) )); then
##						&& (( $(echo "$SNR_Z2 > 3" | bc -l) )) && (( $(echo "$SNR_Z3 > 3" | bc -l) )) 

##						&& (( $(echo "$SNR_N2 > 3" | bc -l) )) && (( $(echo "$SNR_N3 > 3" | bc -l) )) 

##						&& (( $(echo "$SNR_E2 > 3" | bc -l) )) && (( $(echo "$SNR_E3 > 3" | bc -l) )); then
						
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
r ${i} ${file2} ${file3}
hilbert
sqr
write sac ${i/N_/hilsqr_} ${file2/N_/hilsqr_} ${file3/N_/hilsqr_}
r ${i} ${file2} ${file3}
sqr
write sac ${i/N_/sqr_} ${file2/N_/sqr_} ${file3/N_/sqr_}
r ${i/N_/hilsqr_}
addf ${i/N_/sqr_} ${file2/N_/hilsqr_} ${file2/N_/sqr_} ${file3/N_/hilsqr_} ${file3/N_/sqr_}
div 8
mul 2830
div ${delf}
write sac ${i/HZ/env}
rh ${i/HZ/env}
ch lovrok true
ch T0 GMT ${Pyear} ${Pjday} ${Phr} ${Pmin} ${Psec} ${Pmsec}
wh
quit
EOF
## | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g'
var=($(sac << EOF | grep "DEPMEN" | awk '{ print $3 }'
r ${i/HZ/env}
cut T0 -50 -45
r
lh DEPMEN
cut off
r
cut T0 -45 -40
r
lh DEPMEN
cut off
r
cut T0 -40 -35
r
lh DEPMEN
cut off
r
cut T0 -35 -30
r
lh DEPMEN
cut off
r
cut T0 -30 -25
r
lh DEPMEN
cut off
r
cut T0 -25 -20
r
lh DEPMEN
cut off
r
cut T0 -20 -15
r
lh DEPMEN
cut off
r
cut T0 -15 -10
r
lh DEPMEN
cut off
quit
EOF
))
						rm -f ${i/N_/hilsqr_} ${i/N_/sqr_} ${file2/N_/hilsqr_} ${file2/N_/sqr_} ${file3/N_/hilsqr_} ${file3/N_/sqr_}

##						avg1=${var[0]}
##						avg2=${var[1]}
##						avg3=${var[2]}
##						avg4=${var[3]}
##						avg5=${var[4]}
##						avg6=${var[5]}
##						avg7=${var[6]}
##						avg8=${var[7]}
					
						avg_pos=`printf '%s\n' "${var[@]}" |
						awk 'NR==1{min=$1;pos=NR}NR>1 && $1<min{min=$1;pos=NR}END{print pos}'`

##						avg_min=`perl -e "use List::Util qw[min max];print min(${avg1},${avg2},${avg3},${avg4},${avg5},${avg6},${avg7},${avg8})"`
##						echo $avg1 $avg2 $avg3 $avg4 $avg5 $avg6 $avg7 $avg8 $avg_min $avg_pos
##						avg1=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg1"`
##						avg2=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg2"`
##						avg3=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg3"`
##						avg4=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg4"`
##						avg5=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg5"`
##						avg6=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg6"`
##						avg7=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg7"`
##						avg8=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg8"`
##						avg_min=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$avg_min"`

						
						if [ "${avg_pos}" -eq 1 ]; then
##							echo ${avg1}
							diff=$(echo "$absStime - $absPtime + 50" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						elif [ "${avg_pos}" -eq 2 ]; then
##							echo ${avg2}
							diff=$(echo "$absStime - $absPtime + 45" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						elif [ "${avg_pos}" -eq 3 ]; then
##							echo ${avg3}
							diff=$(echo "$absStime - $absPtime + 40" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						elif [ "${avg_pos}" -eq 4 ]; then
##							echo ${avg4}
							diff=$(echo "$absStime - $absPtime + 35" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						elif [ "${avg_pos}" -eq 5 ]; then
##							echo ${avg5}
							diff=$(echo "$absStime - $absPtime + 30" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						elif [ "${avg_pos}" -eq 6 ]; then
##							echo ${avg6}
							diff=$(echo "$absStime - $absPtime + 25" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						elif [ "${avg_pos}" -eq 7 ]; then
## 							echo ${avg7}
							diff=$(echo "$absStime - $absPtime + 20" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						else
##							echo ${avg8}
							diff=$(echo "$absStime - $absPtime + 15" | bc -l)
							strt_ar=($(seq 10 5 100))
							for (( j=0; j<=18; j++ ))
							do
								Sstrt=${strt_ar[${j}]}
								SNR_env=`sac_snr -A-2 -S${Sstrt} -N-${diff} -W5 ${i/HZ/env} | awk '{ print $2 }'`
								if (( $(echo "$SNR_env < 3" | bc -l) )); then
									cut_W=${Sstrt}
									break;
								else
									cut_W=${Sstrt}
								fi
							done
						fi
##						rm -f cut*.ascii
##						rm -f ${i/HZ/Swin} ${i/HZ/Codawin} ${i/HZ/env2}
						temp=${i/HZ/envwin}
sac <<EOF
cut A -1 ${cut_W}
r ${i/HZ/env}
write sac ${i/HZ/envful}
write alpha ${temp/sacii/ascii}
cut off
quit
EOF
##						tail -n +31 ${temp/sacii/ascii} > test.txt
						awk 'NR>30{for (i=1;i<=NF;i++) print $i}' ${temp/sacii/ascii} >> ${h}/outful_N_${stnm}_${comp:0:1}envwin_${frange}.txt
						echo $Sdist >> ${h}/outful_N_${stnm}_${comp:0:1}envwin_${frange}.txt

						rm -f ${temp/sacii/ascii} ${i/HZ/env}
					fi
				fi
			fi
		done
	fi
	rm -f ${h}/N_*HH*.sacii ${h}/N_*Henvful*.ascii ${h}/N_*Henvful*.sacii
##	cd ..
done
#SOD
echo "success"

# for h in /home/user/env_EGEL/husn_events/2011-10-21T11M10S1*
# # /home/user/env_EGEL/husn_events/{2013..2018..1}* /home/user/env_EGEL/husn_events/2019-{01..05..1}*
# do
# 	hyp_f=`find ${h}/ -name "husn*.hyp"`
# 	if [ "${#hyp_f}" -gt 0 ]; then
# 		err_xx=`grep "STATISTICS" ${hyp_f} | awk '{ print $9 }'`
# 		err_yy=`grep "STATISTICS" ${hyp_f} | awk '{ print $15 }'`
# 		err_zz=`grep "STATISTICS" ${hyp_f} | awk '{ print $19 }'`
# 		err_h=$(echo "scale=5;sqrt(${err_xx} + ${err_yy})" | bc -l)
# 		err_zz=$(echo "scale=5;sqrt(${err_zz})" | bc -l)
# 		echo ${err_h} ${err_zz}

# 		if (( $(echo "${err_h} < 5" | bc) )) && (( $(echo "${err_zz} < 5" | bc) )); then

# 			for i in ${h}/N_*envful*.sacii
# 			do
#  				IFS="_"
#  				read -ra f_info <<< "${i}"
#  				IFS=${OLDIFS}
#  				temp1=${f_info[5]}
#  				frange=${temp1//.sacii}
#  ##				echo $frange
#  				stnm=${f_info[3]}
#  				comp=${f_info[4]}

# 		 		O_info=`grep "GEOGRAPHIC" ${hyp_f}`
#  				read -ra O_arr <<< "${O_info}"

# 		 		Odep=${O_arr[13]}

#  				S_info=`grep -w "${stnm}" ${hyp_f} | grep "? S"`
#  				read -ra S_arr <<< "${S_info}"				

# 				Sdist=${S_arr[21]}
# 				if [ "${Sdist}" = "0.0000" ]; then
# 					Sdist="400"
# 				fi
# 				Sdist=$(echo "scale=4;sqrt(${Sdist}^2 + ${Odep}^2)" | bc -l)

# #				Sdist=`tail -n1 ${h}/outful_N_${stnm}_${comp:0:1}envwin_${frange}.txt`				

# 				echo ${i}
# #				echo ${stnm} ${comp} ${frange} ${Sdist}
# sac <<EOF
# r ${i}
# write alpha ${i/sacii/ascii}
# quit
# EOF
# 				awk 'NR>30{for (i=1;i<=NF;i++) print $i}' ${i/sacii/ascii} >> ${h}/outful_N_${stnm}_${comp:0:1}envwin_${frange}.txt
# 				echo $Sdist >> ${h}/outful_N_${stnm}_${comp:0:1}envwin_${frange}.txt
#                 echo "success"

# 			done
# 		else
# 			echo "unwanted files"
# 			rm -f ${h}/N_*envful*.sacii
# 		fi
# 	fi
# done



# #: <<WOD
# matlab <<EOF
# %profile on;
# rw = [dir("/home/user/env_EGEL/husn_events/2014-08-15T04M29S2*");...
#     %dir("/home/user/env_EGEL/husn_events/2012-07-06T08M19S0*");...
#     %dir("/home/user/env_EGEL/husn_events/2012*");...
#     %dir("/home/user/env_EGEL/husn_events/2013*");...
#     %dir("/home/user/env_EGEL/husn_events/2014*");...
#     %dir("/home/user/env_EGEL/husn_events/2015*");...
#     %dir("/home/user/env_EGEL/husn_events/2016*");...
#     %dir("/home/user/env_EGEL/husn_events/2017*");...
#     %dir("/home/user/env_EGEL/husn_events/2018*");...
#     %dir("/home/user/env_EGEL/husn_events/2019*");
#     ];
# mw = struct2cell(rw);
# %length(mw)
# %%{
# for j=1:length(mw)
#     chek = sprintf("%s/%s",mw{(j-1)*6+2},mw{(j-1)*6+1});
#     chek
#     df_1 = sprintf("%s/%s",chek,'outful_N*.txt');
#     alpha = dir(df_1);
#     beta = struct2cell(alpha);
#     %length(alpha)
#     if (~isempty(alpha))
#         for p = 1:length(alpha)
#             if (beta{(p-1)*6+4} > 24700)
#                 df_3 = sprintf("%s/%s",chek,beta{(p-1)*6+1});
#                 %df_3
#                 %A = load(df_3);
#                 fid = fopen(df_3,'r');
#                 temp = textscan(fid,'%f');
#                 fclose(fid);
#                 A = temp{1};
#                 b = bartlett(201);
#                 D = conv(A(1:end-1),b/sum(b),'same');
#                 %pkg load signal
#                 [q,loc] = findpeaks(D(1001:end));
#                 E = 1./D;
#                 [qd,locd] = findpeaks(E(1001:end));
#                 loc = loc + 1000;
#                 locd = locd + 1000;
#                 %length(loc)
#                 %length(locd)
#                 x = 1:1:length(D);
#                 %%plot(x,D); hold on; plot(x(loc),q,'ro*'); plot(x(locd),D(locd),'o'); hold off;
#                 jmp = zeros(length(loc),1);
#                 idx = zeros(length(loc),1);
#                 for k=1:length(loc)
#                     sw = find(locd<loc(k));
#                     jmp(k) = 0;
#                     idx(k) = 0;
#                     if (~isempty(sw))
#                         idx(k) = locd(sw(end));
#                         jmp(k) = D(loc(k))/D(idx(k));
#                     end
#                 end
#                 rw = find(jmp > 4);
#                 cut_p = length(D) - 1;
#                 if (~isempty(rw))
#                     cut_p = idx(rw(1));
#                 end
#                 fid = fopen(df_3,'w');
#                 fprintf(fid,"%e\n",D(1:cut_p)');
#                 fprintf(fid,"%f",A(end));
#                 fclose(fid);
#             end
#         end
#     end
#     %profile off;
#     %T = profile("info");
#     %profshow(T);
# end
# %}
# clearvars;
# quit
# EOF
#WOD
