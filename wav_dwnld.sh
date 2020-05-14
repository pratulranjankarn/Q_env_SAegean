#!/bin/bash

rm -f req_data.txt

OLDIFS=${IFS}

mkdir /home/user/env_EGEL/kef_events

SECONDS=0

# while read V1 V2 V3; do
#     aw=`find /home/user/env_EGEL/HUSN_PZs/* -name "SAC_PZs_*${V1}_HHZ_*"`
#     IFS="_"
#     read -ra nt_arr <<< ${aw}
#     IFS=${OLDIFS}
#     echo ${nt_arr[4]} ${V1} >> stntwrk_husn.txt
# done < stlist_husn.txt

for i in /home/user/NLLOC_locations/loc_kef/husn.2*.hyp
do
    temp=`grep "STATISTICS" ${i}`
    read -ra err_set <<< $temp
    covxx=${err_set[8]}
    covyy=${err_set[14]}
    covzz=${err_set[18]}
    err_hz=$(echo "scale=5;sqrt(${covxx}+${covyy})" | bc -l)
    err_vrt=$(echo "scale=5;sqrt(${covzz})" | bc -l)
# echo ED1
    temp2=`grep "QUALITY" ${i}`
    read -ra rms_set <<< $temp2
    rms=${rms_set[8]}
    if (( $(echo "${err_hz} < 5" | bc) )) && (( $(echo "${err_vrt} < 5" | bc) )); then
        echo $i
        temp=`grep "GEOGRAPHIC" ${i}`
        read -ra evnt_orig <<< $temp
        # echo ${evnt_orig[11]} ${evnt_orig[9]} ${evnt_orig[13]} >> kef_evnts.txt
        yeer=${evnt_orig[2]}
        mnth=${evnt_orig[3]}
        dyt=${evnt_orig[4]}
        hr=${evnt_orig[5]}
        mn=${evnt_orig[6]}
        ss=${evnt_orig[7]}
        IFS="."
        read -ra scmsc <<< $ss
        IFS=${OLDIFS}
        secs=${scmsc[0]}
        if [ "${#secs}" -eq 1 ]; then
            secs=0${secs}
        fi
        mscs=${scmsc[1]}
        echo ${yeer}-${mnth}-${dyt}H${hr}M${mn}S${secs}
        mkdir "/home/user/env_EGEL/kef_events/${yeer}-${mnth}-${dyt}H${hr}M${mn}S${secs}"
        cd "/home/user/env_EGEL/kef_events/${yeer}-${mnth}-${dyt}H${hr}M${mn}S${secs}"

        rm -f req_data.txt
# echo ED2
        day=`date -d ${yeer}/${mnth}/${dyt} +%j`
        startmin=$((10#${mn}-1))
        starthour=$hr
        startday=$day
        startyear=$yeer
        startsec=$secs
        if [ $startmin -lt "0" ]; then
    	    startmin=$((10#${startmin}+60))
    	    starthour=$((10#${starthour}-1))
        fi
        if [ $starthour -lt "0" ]; then
    	    startday=$((10#${day}-1))
    	    starthour=$((10#${starthour}+24))
        fi
        if [ $startday -lt "0" ]; then
    	    startday=$((10#${startday}+365))
    	    startyear=$((10#${yeer}-1))
        fi
        if [ ${#startday} -eq 1 ]; then
    	    startday="0${startday}"
        fi
        if [ ${#startday} -eq 2 ]; then
    	    startday="0${startday}"
        fi
        if [ ${#starthour} -eq 1 ]; then
    	    starthour="0${starthour}"
        fi
        if [ ${#startmin} -eq 1 ]; then
    	    startmin="0${startmin}"
        fi
        if [ ${#secs} -eq 1 ]; then
    	    startsec="0${secs}"
        fi
        jtemp=`echo "${startday} - 1" | bc -l`
        startdm=`date -d "${startyear}/01/01 ${jtemp} days" +%m-%d` 
#        echo $startyear,$startdm,$starthour,$startmin,$startsec
        grep -A100 -P "PHASE" ${i} | grep "S_" | awk '{ print $1 }' > /home/user/env_EGEL/check.txt
# echo ED3
        endmin=$((10#${mn}+4))
        endsec=$secs
        endhour=$hr
        endday=$day
        endyear=$yeer
        if [[ "${endmin#0}" -ge 60 ]]; then
    	    endmin=$((10#${endmin}-60))
    	    endhour=$((10#${endhour}+1))
        fi
        if [[ "${endhour#0}" -ge 24 ]]; then
    	    endday=$((10#${day}+1))
    	    endhour=$((10#${endhour}-24))
        fi
        if [ $endday -gt "365" ]; then
    	    endday=$((10#${endday}-365))
    	    endyear=$((10#${yeer}+1))
        fi
        if [ ${#endday} -eq 1 ]; then
    	    endday="0${endday}"
        fi
        if [ ${#endday} -eq 2 ]; then
    	    endday="0${endday}"
        fi
        if [ ${#endhour} -eq 1 ]; then
    	    endhour="0${endhour}"
        fi
        if [ ${#endmin} -eq 1 ]; then
    	    endmin="0${endmin}"
        fi
        if [ ${#secs} -eq 1 ]; then
        	endsec="0${secs}"
        fi
        jtemp=`echo "${endday} - 1" | bc -l`
        enddm=`date -d "${endyear}/01/01 ${jtemp} days" +%m-%d`
        #echo $endyear,$enddm,$endhour,$endmin,$endsec
        while read V1; do
            nt=`grep -w "${V1}" /home/user/env_EGEL/stntwrk_husn.txt | awk '{ print $1 }'`
            if [ "${#nt}" -gt "0" ]; then
            echo $nt $V1 -- HHZ ${startyear}-${startdm}T${starthour}:${startmin}:${startsec} \
            ${endyear}-${enddm}T${endhour}:${endmin}:${endsec} >> req_data.txt
            echo $nt $V1 -- HHN ${startyear}-${startdm}T${starthour}:${startmin}:${startsec} \
            ${endyear}-${enddm}T${endhour}:${endmin}:${endsec} >> req_data.txt
            echo $nt $V1 -- HHE ${startyear}-${startdm}T${starthour}:${startmin}:${startsec} \
            ${endyear}-${enddm}T${endhour}:${endmin}:${endsec} >> req_data.txt
            fi
        done < /home/user/env_EGEL/check.txt
# echo ED4
#        echo $SECONDS
#         mw=`arclinktool -u Pratul -r req_data.txt eida.gein.noa.gr:80`
# #        echo $SECONDS
#         id=`echo ${mw} | awk '{ print $19 }'`
#         arclinktool -u Pratul -s $id eida.gein.noa.gr:18001 > checking.txt
#         qw=`grep "UNSET" ./checking.txt`
#         while [[ "${#qw}" -gt 0 ]]; do
#             sleep .05
#             echo "repeating"
#             arclinktool -u Pratul -s $id eida.gein.noa.gr:18001 > checking.txt
#             qw=`grep "UNSET" ./checking.txt`
#         done
# #        echo $SECONDS
#         arclinktool -u Pratul -d $id -o data.mseed eida.gein.noa.gr:18001
        curl -v http://eida.gein.noa.gr/fdsnws/dataselect/1/query --data-binary "$(<req_data.txt)" >> data.mseed
        export ALT_RESPONSE_FILE="/home/user/env_EGEL/HUSN_PZs/stations_PZs.dseed"
        rdseed -d -o 1 -E -f data.mseed
        cp $i .
#    rdseed <<EOF
#data.mseed


#d






#Y







#quit
#EOF
       rm -f data.mseed check.txt temp.txt checking.txt
    fi
# echo ED5
#fi
    cd ..
done
echo $SECONDS
