#!/bin/bash

rm -f husn_shlw.pick

for i in /home/user/env_EGEL/husn_events/2010* /home/user/env_EGEL/husn_events/2011-01*
do
	echo $i
	dyt_frst=`find ${i} -name "husn*.pick"`
	IFS='/'
	read -ra arr <<< "${dyt_frst}"
	IFS='.'
	read -ra dyt_scnd <<< "${arr[6]}"
	#echo ${dyt_scnd[1]}
	dyt_fnl=${dyt_scnd[1]}
	tym_o=${dyt_scnd[2]}
	#echo $tym_o
	IFS=' '
	tym_frst=`grep " H=" ${i}/husn*.pick | awk '{ print $2 }'`
	tym_frst=${tym_frst//H=}
	#echo ${tym_frst}
	y_o=${dyt_fnl:0:4}
	m_o=${dyt_fnl:4:2}
	d_o=${dyt_fnl:6:2}
	#echo ${dyt_fnl}
	pick_st=`grep -n -w "STN PHASE" ${i}/husn*.pick | awk '{ print $1 }'`
	pick_st="$(echo -e "${pick_st}" | tr -d '[:space:]')"
	pick_st=${pick_st//:}
	#echo $pick_st
	pick_st=$(echo "${pick_st}+1" | bc)
	pick_nd="-1"
	sed -n "$pick_st,$"p ${i}/husn*.pick > check2.txt
	st_tmp=""
	while read -ra arr 
	do
	#echo ${#arr}
	if [ "${#arr}" -gt 1 ]; then
		#echo ${arr[0]} ${arr[1]} ${arr[2]}
		if [ "${arr[0]:0:2}" = "ES" ]; then
			stnm=${st_tmp}
			phz=${arr[0]}
			IFS=":"
			read -ra tym <<< "${arr[1]}"
			IFS=" "
			if [ "${#tym[0]}" -eq 0 ]; then
				hr_p=${hr_tmp}
				min_p=${min_tmp}
				ss_p=${tym[1]}
			elif [ "${#tym[0]}" -gt 0 ] && [ "${#tym[2]}" -eq 0 ]; then
				hr_p=${hr_tmp}
				min_p=${tym[0]}
				ss_p=${tym[1]}
			else
				hr_p=${tym[0]}
				min_p=${tym[1]}
				ss_p=${tym[2]}
			fi

			if [ "${hr_p:0:1}" -eq 0 ] && [ "${tym_o:0:1}" -eq 2 ]; then
				d_o=$(echo "${d_o}+1" | bc)
			fi
			m_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
			n_d=`cal ${m_o} ${y_o} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
			if [ "${#m_chk}" -eq 0 ]; then
				m_o=$(echo "${m_o} + 1" | bc)
				d_o=$(echo "${d_o}-${n_d}" | bc)
			fi
			y_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
			if [ "${#y_chk}" -eq 0 ]; then
				y_o=$(echo "${y_o} + 1" | bc)
				m_o=$(echo "${m_o} - 12" | bc)
			fi
			dyt_fnl=${y_o}${m_o}${d_o} 

			#echo ${stnm} ${arr[0]} ${hr_p} ${min_p} ${ss_p}
		else
			st_tmp=${arr[0]}
			stnm=${st_tmp}
			phz=${arr[1]}
			IFS=":"
			read -ra tym <<< "${arr[2]}"
			IFS=" "
			hr_p=${tym[0]}
			min_p=${tym[1]}
			ss_p=${tym[2]}
			hr_tmp=${hr_p}
			min_tmp=${min_p}
			#echo ${hr_p} ${tym_o}
			if [ "${hr_p:0:1}" -eq 0 ] && [ "${tym_o:0:1}" -eq 2 ]; then
				d_o=$(echo "${d_o}+1" | bc)
			fi
			m_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
			n_d=`cal ${m_o} ${y_o} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
			if [ "${#m_chk}" -eq 0 ]; then
				m_o=$(echo "${m_o} + 1" | bc)
				d_o=$(echo "${d_o}-${n_d}" | bc)
			fi
			y_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
			if [ "${#y_chk}" -eq 0 ]; then
				y_o=$(echo "${y_o} + 1" | bc)
				m_o=$(echo "${m_o} - 12" | bc)
			fi
			dyt_fnl=${y_o}${m_o}${d_o} 

			#echo ${arr[0]} ${arr[1]} ${hr_p} ${min_p} ${ss_p}
		fi

	hrm_fnl=${hr_p}${min_p}
	ss_fnl=${ss_p}00
	if [ "${phz:2:1}" = "G" ]; then
		ph_err="0.05"
		phz=${phz:1:1}_0
	elif [ "${phz:2:1}" = "B" ]; then
		ph_err="0.10"
		phz=${phz:1:1}_1
	else
		ph_err="0.15"
		phz=${phz:1:1}_2
	fi
	# ph_err=`printf "%0.2f\n" ${ph_err}`
	echo ${stnm} ? E ? ${phz} ? ${dyt_fnl} ${hrm_fnl} ${ss_fnl} GAU ${ph_err}e+00 -1.00e+00 -1.00e+00 -1.00e+00 >> husn_shlw.pick
	fi
	done < check2.txt
	echo " " >> husn_shlw.pick
done


for i in /home/user/env_EGEL/husn_events/2011-{02..12..1}*  /home/user/env_EGEL/husn_events/{2012..2018..1}*  /home/user/env_EGEL/husn_events/2019-01*  /home/user/env_EGEL/husn_events/2019-02-{01..04..1}*  
do
echo $i
dyt_frst=`grep -w -A1 "Date" ${i}/husn*.pick | tail -n1 | awk '{ print $1 }'`
tym_frst=`grep -w -A1 "Date" ${i}/husn*.pick | tail -n1 | awk '{ print $2 }'`
y_o=${dyt_frst:0:4}
m_o=${dyt_frst:5:2}
d_o=${dyt_frst:8:2}
#echo ${dyt_fnl}
pick_st=`grep -n -w "Sta     Dist" ${i}/husn*.pick | awk '{ print $1 }'`
pick_st="$(echo -e "${pick_st}" | tr -d '[:space:]')"
pick_st=${pick_st//:Sta}
pick_st=$(echo "${pick_st}+1" | bc)
pick_nd=`grep -n -w "STOP" ${i}/husn*.pick | awk '{ print $1 }'`
pick_nd="$(echo -e "${pick_nd}" | tr -d '[:space:]')"
pick_nd=${pick_nd//:STOP}
if [ "${#pick_nd}" -eq 0 ]; then
sed -n "$pick_st,$"p ${i}/husn*.pick > check.txt
else
pick_nd=$(echo "${pick_nd}-3" | bc)
sed -n "$pick_st,$pick_nd"p ${i}/husn*.pick > check.txt
fi
while read stnm ep_dst azm phz hrminss t_res T_T qual comp
do
if [ "${phz:0:1}" = "P" ] || [ "${phz:0:1}" = "S" ]; then
stnm="$(echo -e "${stnm}" | tr -d '[:space:]')"
if [ "${hrminss:0:1}" -eq 0 ] && [ "${tym_frst:0:1}" -eq 2 ]; then
d_o=$(echo "${d_o}+1" | bc)
echo ${y_o}${m_o}${d_o}
fi
m_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
n_d=`cal ${m_o} ${y_o} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
if [ "${#m_chk}" -eq 0 ]; then
m_o=$(echo "${m_o} + 1" | bc)
d_o=$(echo "${d_o}-${n_d}" | bc)
echo ${y_o}${m_o}${d_o}
fi
y_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
if [ "${#y_chk}" -eq 0 ]; then
y_o=$(echo "${y_o} + 1" | bc)
m_o=$(echo "${m_o} - 12" | bc)
echo ${y_o}${m_o}${d_o}
fi
dyt_fnl=${y_o}${m_o}${d_o}
hrm_fnl=${hrminss:0:2}${hrminss:3:2}
ss_fnl=${hrminss:6:5}00
if [ ${#phz} -gt 1 ]; then
ph_err=$(echo "scale=2;0.05 + ${phz:2:1}*0.05" | bc -l)
ph_err=`printf "%0.2f\n" ${ph_err}`
else
ph_err="0.05"
fi
echo ${stnm} ? ${comp:2:1} ? ${phz} ? ${dyt_fnl} ${hrm_fnl} ${ss_fnl} GAU ${ph_err}e+00 -1.00e+00 -1.00e+00 -1.00e+00 >> husn_shlw.pick
fi
done < check.txt
echo " " >> husn_shlw.pick
done



for i in /home/user/env_EGEL/husn_events/2019-02-{05..28..1}*  /home/user/env_EGEL/husn_events/2019-{03..05..1}*  
do
echo $i
dyt_frst=`grep -A1 "Date" ${i}/husn*.pick | tail -n1 | awk '{ print $1 }'`
tym_frst=`grep -A1 "Date" ${i}/husn*.pick | tail -n1 | awk '{ print $2 }'`
y_o=${dyt_frst:0:4}
m_o=${dyt_frst:5:2}
d_o=${dyt_frst:8:2}
echo ${dyt_fnl}
pick_st=`grep -n " Sta     Dist" ${i}/husn*.pick | awk '{ print $1 }'`
pick_st="$(echo -e "${pick_st}" | tr -d '[:space:]')"
pick_st=${pick_st//:}
pick_st=$(echo "${pick_st}+1" | bc)
pick_nd=`grep -n -w "STOP" ${i}/husn*.pick | awk '{ print $1 }'`
pick_nd="$(echo -e "${pick_nd}" | tr -d '[:space:]')"
pick_nd=${pick_nd//:}
echo ${pick_st} ${pick_nd}
if [ "${#pick_nd}" -eq 0 ]; then
sed -n "$pick_st,$"p ${i}/husn*.pick > check2.txt
else
pick_nd=$(echo "${pick_nd}-2" | bc)
sed -n "$pick_st,$pick_nd"p ${i}/husn*.pick > check2.txt
fi
while read -ra arr #stnm ep_dst azm phz hrminss t_res T_T qual comp
do
stnm=${arr[0]}
phz=${arr[3]}
hrminss=${arr[4]}
comp=${arr[11]}
if [ "${#comp}" -eq 0 ]; then
	comp=${arr[8]}
fi
echo ${comp}
if [ "${phz:0:1}" = "P" ] || [ "${phz:0:1}" = "S" ]; then
stnm="$(echo -e "${stnm}" | tr -d '[:space:]')"
if [ "${hrminss:0:1}" -eq 0 ] && [ "${tym_frst:0:1}" -eq 2 ]; then
d_o=$(echo "${d_o}+1" | bc)
echo ${y_o}${m_o}${d_o}
fi
m_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
n_d=`cal ${m_o} ${y_o} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
if [ "${#m_chk}" -eq 0 ]; then
m_o=$(echo "${m_o} + 1" | bc)
d_o=$(echo "${d_o}-${n_d}" | bc)
echo ${y_o}${m_o}${d_o}
fi
y_chk=`date -d ${y_o}-${m_o}-${d_o} +%m`
if [ "${#y_chk}" -eq 0 ]; then
y_o=$(echo "${y_o} + 1" | bc)
m_o=$(echo "${m_o} - 12" | bc)
echo ${y_o}${m_o}${d_o}
fi
dyt_fnl=${y_o}${m_o}${d_o}
hrm_fnl=${hrminss:0:2}${hrminss:3:2}
ss_fnl=${hrminss:6:5}00
if [ ${#phz} -gt 1 ]; then
ph_err=$(echo "scale=2;0.05 + ${phz:2:1}*0.05" | bc -l)
ph_err=`printf "%0.2f\n" ${ph_err}`
else
ph_err="0.05"
fi
echo ${stnm} ? ${comp:5:1} ? ${phz} ? ${dyt_fnl} ${hrm_fnl} ${ss_fnl} GAU ${ph_err}e+00 -1.00e+00 -1.00e+00 -1.00e+00 >> husn_shlw.pick
fi
done < check2.txt
echo " " >> husn_shlw.pick
done
