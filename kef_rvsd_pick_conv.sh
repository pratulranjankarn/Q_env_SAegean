#!/bin/bash

OLDIFS=${IFS}

rm -f kef_shlw.pick

while read -ra Ev_arr
do
#   echo ${Ev_arr[0]}
    IFS="/"
    read -ra dte <<< ${Ev_arr[0]}
#   echo ${dte[2]}
    IFS=${OLDIFS}

    sed -n "/${dte[0]}\/${dte[1]}\/${dte[2]} ${Ev_arr[1]}/,/DATA_TYPE BULLETIN/p" \
    kef_pick_fnl.txt | grep "T__" > test.txt
    echo ${Ev_arr[0]} ${Ev_arr[1]}
    while read -ra pick_arr 
    do
    y_o=${dte[0]}
    m_o=${dte[1]}
    d_o=${dte[2]}

        if [ "${pick_arr[3]:0:1}" == "P" ] || [ "${pick_arr[3]:0:1}" == "S" ]; then
            stnm=${pick_arr[0]}
            comp=${pick_arr[8]:2:1}
            phz=${pick_arr[3]}
            if [ "${phz:2:1}" = "0" ]; then
		        ph_err="0.05"
	        elif [ "${phz:2:1}" = "1" ]; then
		        ph_err="0.10"
	        elif [ "${phz:2:1}" = "2" ]; then
		        ph_err="0.15"
            else
                ph_err="0.20"
	        fi
            hrminss=${pick_arr[4]}
            tym_frst=${Ev_arr[1]}
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
            echo ${stnm} ? ${comp} ? ${phz} ? ${dyt_fnl} ${hrm_fnl} ${ss_fnl} GAU ${ph_err}e+00 -1.00e+00 -1.00e+00 -1.00e+00 >> kef_shlw.pick
        fi
    done < test.txt
    echo " " >> kef_shlw.pick
done < kef_evinfo_fnl.txt