#!/bin/bash

#rm -rf husn_events

#mkdir husn_events
: <<SOD
while read V1 V2 V3 V4 V5 V6 V7 V8;
do
mkdir ~/env_EGEL/husn_events/${V3}T${V4}
year=${V3:0:4}
mnth=${V3:5:2}
dayt=${V3:8:2}
hwr=${V4:0:2}
min=${V4:3:2}
sec=${V4:6:2}
wget -O- http://bbnet.gein.noa.gr/alerts_manual/${year}/${mnth}/${V2}_info.html | sed '/<*>/d' > ~/env_EGEL/husn_events/${V3}T${V4}/husn.${year}${mnth}${dayt}.${hwr}${min}${sec}.pick
sed -i '1,21d' ~/env_EGEL/husn_events/${V3}T${V4}/husn.${year}${mnth}${dayt}.${hwr}${min}${sec}.pick
done < husn_catalog7
SOD

#: <<EOF
while read V1 V2 V3 V4 V5 V6 V7 V8;
do
mkdir ~/env_EGEL/husn_events/${V3}T${V4}
year=${V3:0:4}
mnth=${V3:5:2}
dayt=${V3:8:2}
hwr=${V4:0:2}
min=${V4:3:2}
sec=${V4:6:2}
if [[ "${V2:0:3}" == "noa" ]]; then
wget -O- http://bbnet.gein.noa.gr/Events/${year}/${mnth}/${V2}_info.html | sed '/<*>/d' > ~/env_EGEL/husn_events/${V3}T${V4}/husn.${year}${mnth}${dayt}.${hwr}${min}${sec}.pick
else
wget -O- http://bbnet.gein.noa.gr/alerts_manual/${year}/${mnth}/${V2}_info.html | sed '/<*>/d' > ~/env_EGEL/husn_events/${V3}T${V4}/husn.${year}${mnth}${dayt}.${hwr}${min}${sec}.pick
fi
done < kef_catalog_2010.txt
#EOF

