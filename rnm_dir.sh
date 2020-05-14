#!/bin/bash

for i in /home/user/env_EGEL/husn_events/2*
do
IFS=':'
read -ra arr <<< ${i}
#echo ${arr[0]}\:${arr[1]}\:${arr[2]} ${arr[0]}M${arr[1]}S${arr[2]}
if [ -a ${arr[0]}\:${arr[1]}\:${arr[2]} ]; then
	mv -v ${arr[0]}\:${arr[1]}\:${arr[2]} ${arr[0]}M${arr[1]}S${arr[2]}
fi
IFS=' '
done
