#!/bin/bash

cd temp_extract2

find . -size -65k -exec rm -f {} \;

for i in 2*
do
ystrt=${i:0:4}
jdstrt=${i:5:3}
jtemp=`echo "${jdstrt} - 1" | bc -l`
mdstrt=`date -d "${ystrt}-01-01 ${jtemp} days" +%m/%d`
hrstrt=${i:9:2}
minstrt=${i:12:2}
scstrt=${i:15:2}
mscstrt=${i:18:3}

#echo ${ystrt}/${mdstrt}T${hrstrt}:${minstrt}:${scstrt}.${mscstrt}

sabstrt=`date -d ${ystrt}/${mdstrt}T${hrstrt}:${minstrt}:${scstrt}.${mscstrt} +%s`

yend=${i:23:4}
jdend=${i:28:3}
jtemp=`echo "${jdend} - 1" | bc -l`
mdend=`date -d "${yend}-01-01 ${jtemp} days" +%m/%d`
mdend2=`date -d "${yend}-01-01 ${jtemp} days" +%m-%d`
#echo $mdend3
hrend=${i:32:2}
minend=${i:35:2}
scend=${i:38:2}
mscend=${i:41:3}

#echo ${yend}/${mdend}T${hrend}:${minend}:${scend}.${mscend}

sabsend=`date -d ${yend}/${mdend}T${hrend}:${minend}:${scend}.${mscend} +%s`

aw=`find ~/env_EGEL/husn_events/2* -name "${yend}-${mdend2}*"`

for j in $aw
do
#echo $j
yrf=${j:32:4}
mnf=${j:37:2}
dyf=${j:40:2}
hrf=${j:43:2}
minf=${j:46:2}
scf=${j:49:2}
#echo ${yrf}/${mnf}/${dyf}T${hrf}:${minf}:${scf}

sabsfdr=`date -d ${yrf}/${mnf}/${dyf}T${hrf}:${minf}:${scf} +%s`

#echo $sabstrt $sabsfdr $sabsend

if [[ "${sabsfdr}" -ge "${sabstrt}" ]] && [[ "${sabsfdr}" -le "${sabsend}" ]]; then
echo $i
echo $j
mv $i $j/
break;
fi
done

done

cd ..
