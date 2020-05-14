#!/bin/bash

for i in /home/user/env_EGEL/husn_events/2019-05-27T12:50:18
do
echo $i
cd "${i}"
find ./ -name "out*envwin_2-4.sacii" > /home/user/env_EGEL/test.txt
if [ -s /home/user/env_EGEL/test.txt ]; then
while read V1
do


done < /home/user/env_EGEL/test.txt
fi
cd ../..
done
