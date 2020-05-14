#!/bin/bash

awk '{ print $3"T"$4 }' husn_catalog5 > select_files.txt

while read i
do
echo $i 
cd "./husn_events/${i}"
rm -f req_data.txt
echo "REQUEST WAVEFORM format=MSEED" >> req_data.txt
year=${i:0:4}
mnth=${i:5:2}
dte=${i:8:2}
hwrchek=${i:11:2}
hwr=${i:11:2}
min=${i:14:2}
sec=${i:17:2}
#echo ${i:14:4} ${i:19:2} ${i:22:2} ${i:25:2} ${i:28:2} ${i:31:2}
day=`date -d ${year}/${mnth}/${dte} +%j`
startmin=$((10#${min}-1))
starthour=$hwr
startday=$day
startyear=$year
startsec=$sec
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
    	startyear=$((10#${year}-1))
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
    if [ ${#sec} -eq 1 ]; then
    	startsec="0${sec}"
    fi
jtemp=`echo "${startday} - 1" | bc -l`
startdm=`date -d "${startyear}/01/01 ${jtemp} days" +%m,%d` 
head -n100 ./husn.20* | grep " EP" > check.txt
while read V1 V2 V3 V4 V5 V6 V7 V8 V9;
do
Y1=`grep -A1 "${V1}" ./husn.20* | grep " ES"`
if [[ "${#Y1}" -gt 0 ]]; then
hwr=${V3:0:2}
min=${V3:3:2}
sec=${V3:6:2}
endmin=$((10#${min}+3))
endsec=$sec
endhour=$hwr
endday=$day
if [[ "${hwr}" -eq "00" ]] && [[ "${hwrchek}" -ne "00" ]]; then
endday=$((10#${day}+1))
fi
endyear=$year
    if [[ "${endmin#0}" -ge 60 ]]; then
    	endmin=$((10#${endmin}-60))
    	endhour=$((10#${endhour}+1))
    fi
    if [[ "${endhour#0}" -ge 24 ]]; then
    	endday=$((10#${day}+1))
    	endhour=$((10#${endhour}-24))
    fi
    if [ ${endday} -gt "365" ]; then
    	endday=$((10#${endday}-365))
    	endyear=$((10#${year}+1))
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
    if [ ${#sec} -eq 1 ]; then
    	endsec="0${sec}"
    fi
jtemp=`echo "${endday} - 1" | bc -l`
enddm=`date -d "${endyear}/01/01 ${jtemp} days" +%m,%d`
echo ${startyear},${startdm},${starthour},${startmin},${startsec} ${endyear},${enddm},${endhour},${endmin},${endsec} _HUSN ${V1} HHZ . >> req_data.txt
echo ${startyear},${startdm},${starthour},${startmin},${startsec} ${endyear},${enddm},${endhour},${endmin},${endsec} _HUSN ${V1} HHN . >> req_data.txt
echo ${startyear},${startdm},${starthour},${startmin},${startsec} ${endyear},${enddm},${endhour},${endmin},${endsec} _HUSN ${V1} HHE . >> req_data.txt
fi

done < check.txt
echo "END" >> req_data.txt
#echo $SECONDS
mw=`arclinktool -u Pratul -r req_data.txt eida.gein.noa.gr:18001`
#echo $SECONDS
id=`echo ${mw} | awk '{ print $19 }'`
arclinktool -u Pratul -s $id eida.gein.noa.gr:18001 > checking.txt
qw=`grep "UNSET" ./checking.txt`
while [[ "${#qw}" -gt 0 ]]; do
sleep .05
echo "repeating"
arclinktool -u Pratul -s $id eida.gein.noa.gr:18001 > checking.txt
qw=`grep "UNSET" ./checking.txt`
done
#echo $SECONDS
arclinktool -u Pratul -d $id -o data.mseed eida.gein.noa.gr:18001
export ALT_RESPONSE_FILE="/home/user/env_EGEL/HUSN_PZs/Package_1557982854801-NOA_832811_dataless.dseed"
rdseed -d -o 1 -E -f data.mseed
rm -f data.mseed check.txt temp.txt checking.txt
cd ../..
done < select_files.txt
