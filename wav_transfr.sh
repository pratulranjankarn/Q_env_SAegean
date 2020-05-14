#!/bin/bash

mkdir /run/media/user/08B7-FB2E/EGELADOS_data/

for i in ./EGEL_events/*
do
echo $i
mkdir /run/media/user/08B7-FB2E/EGELADOS_data/${i:14:19}
cp -r ${i}/*.SAC ${i}/*.hyp /run/media/user/08B7-FB2E/EGELADOS_data/${i:14:19}/

done
