#!/bin/bash

OLDIFS=${IFS}

for h in /home/user/env_EGEL/EGEL_events/20*
#/home/user/env_EGEL/EGEL_events/2010* /home/user/env_EGEL/EGEL_events/2011-{01..07..1}* 
do
	echo "${h}"
	#sudo swapoff -a
	#sudo swapon -a	
	#pwd
	rm -f ${h}/*.sacii
	hyp_f=`find ${h}/ -name "egel*.hyp"`
	if [ "${#hyp_f}" -gt 0 ]; then
		err_xx=`grep "STATISTICS" ${hyp_f} | awk '{ print $9 }'`
		err_yy=`grep "STATISTICS" ${hyp_f} | awk '{ print $15 }'`
		err_zz=`grep "STATISTICS" ${hyp_f} | awk '{ print $19 }'`
		err_h=$(echo "scale=5;sqrt(${err_xx} + ${err_yy})" | bc -l)
		err_zz=$(echo "scale=5;sqrt(${err_zz})" | bc -l)
		echo ${err_h} ${err_zz}

		if (( $(echo "${err_h} < 5" | bc) )) && (( $(echo "${err_zz} < 5" | bc) )); then
			for i in ${h}/2*.*H*
			do
#				echo $i
				if [[ "${#i}" -gt 64 ]]; then
					chek_sz=`wc -c ${i} | awk '{ print $1 }'`
					if [ "${chek_sz}" -gt 70000 ] && [ "${chek_sz}" -lt 170000 ]; then
						echo $i
						IFS="/"
						read -ra adr <<< "${i}"
						IFS="."
						read -ra arr <<< "${adr[6]}"
						stnm=${arr[2]}
						comp=${arr[3]}
#						echo $stnm $comp
						IFS=${OLDIFS}
						aw=`grep -w "$stnm" ${hyp_f} | grep " S"`
	
						if [[ "${#aw}" -gt 0 ]]; then
							wav_year=${arr[0]:0:4}
							wav_jday=${arr[0]:4:3}
							jtemp=`echo "${wav_jday} - 1" | bc -l`
							wav_md=`date -d "${wav_year}-01-01 ${jtemp} days" +%m/%d`
							wav_hwr=${arr[0]:7:2}
							wav_min=${arr[0]:9:2}
							wav_sec=${arr[0]:11:2}
							wav_msec=${arr[1]}
#							echo ${wav_year}/${wav_md}T${wav_hwr}:${wav_min}:${wav_sec}.${wav_msec}
							pzf=`find /home/user/env_EGEL/EGEL_PZs/* -name "SAC*${stnm}_${comp}*"`
#							echo $pzf 
							j=0
							for j in ${pzf};
							do
								#echo $j
								IFS="_"
								read -ra pzar <<< "${j}"
								#echo ${pzar1[0]} ${pzar1[1]} ${pzar1[2]} ${pzar1[3]} ${pzar1[4]} ${pzar1[5]} ${pzar1[6]} ${pzar1[7]} ${pzar1[8]}
								IFS="."
								read -ra pzst <<< "${pzar[8]}"
								read -ra pznd <<< "${pzar[9]}"
								IFS=${OLDIFS}
								A_j=${pznd[1]}
								A_y=${pznd[0]}
								A_h=${pznd[2]}
								A_m=${pznd[3]}
								A_s=${pznd[4]}
								A_h=$(echo "$A_h-1" | bc)
								A_m=$(echo "$A_m-1" | bc)
								A_s=$(echo "$A_s-1" | bc)
								A_us=${pznd[5]}
								A_md=`date -d "${A_y}-01-01 ${A_j} days" +%m/%d`
#								echo ${A_y}/${A_md}T${A_h}:${A_m}:${A_s}
								A_pzf=`date -d "${A_y}/${A_md}T${A_h}:${A_m}:${A_s}" +%s`
		
								B_j=${pzst[1]}
								B_y=${pzst[0]}
								B_h=${pzst[2]}
								B_m=${pzst[3]}
								B_s=${pzst[4]}
								B_us=${pzst[5]}
								B_j=$(echo "$B_j-1" | bc)
								B_md=`date -d "${B_y}-01-01 ${B_j} days" +%m/%d`
								B_pzf=`date -d "${B_y}/${B_md}T${B_h}:${B_m}:${B_s}.${B_us}" +%s`
#								echo ${B_y}/${B_md}T${B_h}:${B_m}:${B_s}.${B_us}
		
								A_o=`date -d "${wav_year}/${wav_md}T${wav_hwr}:${wav_min}:${wav_sec}.${wav_msec}" +%s`
								#echo $A_o $A_pzf $B_pzf
								if [[ "${A_o}" -ge "${B_pzf}" ]] && [[ "${A_o}" -le "${A_pzf}" ]]; then
									o_var=$j
								fi
							done
							pzf=$o_var
#							echo $pzf
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
		
#							rm -f ${h}/N_${stnm}_${comp}_1-2.sacii ${h}/${stnm}_${comp}_1-2.ascii
#							rm -f ${h}/N_${stnm}_${comp}_2-4.sacii ${h}/${stnm}_${comp}_2-4.ascii
#							rm -f ${h}/N_${stnm}_${comp}_4-8.sacii ${h}/${stnm}_${comp}_4-8.ascii
#							rm -f ${h}/N_${stnm}_${comp}_8-16.sacii ${h}/${stnm}_${comp}_8-16.ascii
#							rm -f ${h}/N_${stnm}_${comp}_16-32.sacii ${h}/${stnm}_${comp}_16-32.ascii
#							rm -f ${h}/N_${stnm}_${comp}.sacii ${h}/${stnm}_${comp}.ascii
		
		
							#echo "S1"
sac <<EOF
wild echo off
read ${i}
ch lovrok true
wh
rmean
rtrend
trans from polezero s ${pzf} to vel freqlimits ${c1} ${c2} ${c3} ${c4}
write SAC ${h}/N_${stnm}_${comp}.sacii
dc ALL
quit
EOF


							if [[ "${comp:0:1}" = "B" ]]; then
sac <<KOD
wild echo off
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 1 2
w SAC ${h}/N_${stnm}_${comp}_1-2.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 2 4
w SAC ${h}/N_${stnm}_${comp}_2-4.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 4 8
w SAC ${h}/N_${stnm}_${comp}_4-8.sacii
dc ALL
quit
KOD
							elif [[ "${comp:0:1}" = "S" ]]; then
sac <<MOD
wild echo off
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 1 2
w SAC ${h}/N_${stnm}_${comp}_1-2.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 2 4
w SAC ${h}/N_${stnm}_${comp}_2-4.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 4 8
w SAC ${h}/N_${stnm}_${comp}_4-8.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 8 16
w SAC ${h}/N_${stnm}_${comp}_8-16.sacii
dc ALL
quit
MOD
							else
#echo "S2"
sac <<SOD
wild echo off
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 1 2
w SAC ${h}/N_${stnm}_${comp}_1-2.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 2 4
w SAC ${h}/N_${stnm}_${comp}_2-4.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 4 8
w SAC ${h}/N_${stnm}_${comp}_4-8.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 8 16
w SAC ${h}/N_${stnm}_${comp}_8-16.sacii
r ${h}/N_${stnm}_${comp}.sacii
bp bu p 2 n 2 c 16 32
w SAC ${h}/N_${stnm}_${comp}_16-32.sacii
dc ALL
quit
SOD
							fi
							rm -f ${h}/N_${stnm}_${comp}.sacii ${h}/N_${stnm}_${comp}.ascii
							#: << HOLD
							#HOLD
						fi
	
					fi
				fi
			done
		fi
	fi
done
