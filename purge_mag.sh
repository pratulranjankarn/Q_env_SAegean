grep -E "ML     3.6|ML     3.7|ML     3.8|ML     3.9|ML     4.*" kef_evinfo_fnl.txt\
 | awk '{ print $1"H"$2 }' | sed 's\/\-\g;s/\..*//g;s/:/M/1;s/:..//g' > test.txt


rm -f test2.txt

while read V1; do
# echo $V1
V2=`find ./kef_events/ -name "*$V1*" -type d`
if [ "${#V2}" -gt 0 ]; then
rm -rf $V2
echo $V2
fi
done < test.txt

grep -A7 "DATA_TYPE" cornoth_pick_fnl.txt \
| grep -E -w -B3 "ML     3.6*|ML     3.7*|ML     3.8*|ML     3.9*|ML     4.*" \
 | grep "m i ke NOA-GEIN" | awk '{ print $1"H"$2 }' | sed 's\/\-\g;s/\..*//g;s/:/M/1;s/:..//g' > test.txt


rm -f test2.txt

while read V1; do
# echo $V1
V2=`find ./cornoth_events/ -name "*$V1*" -type d`
if [ "${#V2}" -gt 0 ]; then
rm -rf $V2
echo $V2
fi
done < test.txt