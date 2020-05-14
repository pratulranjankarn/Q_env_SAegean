#!/bin/bash

stnm="ATH"

rm -f R_list.txt

find ./husn_events -name "env_param.txt" -exec grep -w "ATH" {} + > check.txt

while read V1 V2 V3; 
do
    IFS="/"
    read -ra arr <<< ${V1}
    IFS=" "
    aw=`cat ${arr[0]}/${arr[1]}/${arr[2]}/env_param.txt | wc -l`

    if [ "${aw}" -gt 1 ]; then
        bse=`grep -w "ATH" ${arr[0]}/${arr[1]}/${arr[2]}/env_param.txt | awk '{ print $5}'`
        bse=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$bse"`
        while read A1 A2 A3 A4 A5 A6; do
            if [ "${A1}" != "ATH" ]; then
                A5=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$A5"`
                R_A=`echo "${A5}/${bse}" | bc -l`
                stltln=`grep -w "${A1}" stlist_f.txt | awk '{ print $2,$3 }'`
                echo ${A1} ${stltln} ${R_A} >> R_list.txt
            fi
        done < ${arr[0]}/${arr[1]}/${arr[2]}/env_param.txt
    fi

done < check.txt
