#!/bin/bash

rm -f catalogue_final2.csv

while read V1 V2 V3; do a1=`grep "Lat ${V2} Long ${V1} Depth ${V3}" *.hyp`; read -ra a2 <<<${a1};
    if [ "${a2[7]:1:1}" == "." ]; then
        sec=${a2[7]:0:1}0
    else
        sec=${a2[7]:0:2}
    fi
    echo ${a2[2]}-${a2[3]}-${a2[4]}H${a2[5]}M${a2[6]}S${sec},${V2},${V1},${V3} >> catalogue_final2.csv
    echo ${a2[2]}-${a2[3]}-${a2[4]}H${a2[5]}M${a2[6]}S${sec}
done < seismicity_final.txt