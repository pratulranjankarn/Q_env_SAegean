#!/bin/bash

#: <<SOD
# rm -f Herr.png Herr.txt Verr.txt Verr.png Rmser.png Rmser.txt Deptd.txt Deptd.png Terr.txt Depxy.txt Depxy.png recheck.txt Depxv.png Depxh.png Herr.ps Verr.ps Rmser.ps error_husn.ps depth_dist_husn.ps

# for i in `find ./husn_events/ -mindepth 2 -maxdepth 2 -name "*.hyp"` `find ./EGEL_events/ -mindepth 2 -maxdepth 2 -name "*.hyp"` \
# `find ./Cycnet_events/ -mindepth 2 -maxdepth 2 -name "*.hyp"` `find ./kef_events/ -mindepth 2 -maxdepth 2 -name "*.hyp"` \
# `find ./cornoth_events/ -mindepth 2 -maxdepth 2 -name "*.hyp"`
# do
# depchek=`grep "GEOGRAPHIC*" $i | awk '{print $14}'`
# if (( $(echo "${depchek} <= 40" | bc -l) )); then
# echo $depchek
# CovXX=`grep "CovXX*" $i | awk '{print $9}'`
# CovYY=`grep "CovXX*" $i | awk '{print $15}'`
# CovZZ=`grep "CovXX*" $i | awk '{print $19}'`
# Herr=$(echo "sqrt(${CovXX}+${CovYY})" | bc -l)
# Verr=$(echo "sqrt(${CovZZ})" | bc -l)
# echo $Herr >> Herr.txt
# echo $Verr >> Verr.txt
# Rmser=`grep "QUALITY*" $i | awk '{print $9}'`
# echo $Rmser >> Rmser.txt
# echo $depchek >> Deptd.txt
# fi
# done
#SOD
#paste Deptd.txt Verr.txt > Depxv.txt
#paste Deptd.txt Herr.txt > Depxh.txt

PS1="error_evnts.ps"
PS2="depth_dist_evnts.ps"

rm -f $PS1 $PS2

#: <<EOF

gmt set PS_PAGE_ORIENTATION = Landscape

gmt set MAP_FRAME_PEN = thickest,black

gmt set FONT_LABEL = 14p,Helvetica-Bold,black

gmt set FONT_ANNOT_PRIMARY = 14p,Helvetica-Bold,black

herr_avg=`awk '{if ($1<20) { x+=$1; y+=1;} next} END { printf("%.2f",x/y)}' Herr.txt`
herr_std=`awk '{if ($1<20) { x+=($1-var)^2; y+=1;} next} END { printf("%.2f",sqrt(x/(y-1)))}' var=$herr_avg Herr.txt`

verr_avg=`awk '{if ($1<20) { x+=$1; y+=1;} next} END { printf("%.2f",x/y)}' Verr.txt`
verr_std=`awk '{if ($1<20) { x+=($1-var)^2; y+=1;} next} END { printf("%.2f",sqrt(x/(y-1)))}' var=$verr_avg Verr.txt`

rmser_avg=`awk '{if ($1<1) { x+=$1; y+=1;} next} END { printf("%.2f",x/y)}' Rmser.txt`
rmser_std=`awk '{if ($1<1) { x+=($1-var)^2; y+=1;} next} END { printf("%.2f",sqrt(x/(y-1)))}' var=$rmser_avg Rmser.txt`

gmt pshistogram Herr.txt -JX4.5i/3i -T0/20/1 -W0.5p -Bxf0.5a2+l"Horizontal Error (in Km)" -By+l"Counts" -BWSne -G0/128/128 -K -V -X1.0i -Y4.8i >> $PS1
# -R0/50/0/1200
echo "Mean = ${herr_avg} km" | gmt pstext -R -J -DJ1c/1c -F+cTR -O -K -V >> $PS1
echo "Std = ${herr_std} km" | gmt pstext -R -J -DJ-6.2c/2c -F+cTR+jLT -O -K -V >> $PS1
gmt pshistogram Verr.txt -JX4.5i/3i -T0/20/1 -W0.5p -Bxf0.5a2+l"Depth Error (in Km)" -By+l"Counts" -BWSne -G0/128/128 -K -O -V -X5.8i -Y0.0i >> $PS1
# -R0/50/0/800
echo "Mean = ${verr_avg} km" | gmt pstext -R -J -DJ1c/1c -F+cTR -O -K -V >> $PS1
echo "Std = ${verr_std} km" | gmt pstext -R -J -DJ-6.1c/2c -F+cTR+jLT -O -K -V >> $PS1
gmt pshistogram Rmser.txt -JX4.5i/3i -T0/1/0.05 -W0.5p -Bxf0.02a0.1+l"RMS Residuals (in s)" -By+l"Counts" -BWSne -G0/128/128 -O -V -K -X-5.8i -Y-4.0i >> $PS1
# -R0/1/0/150
echo "Mean = ${rmser_avg} s" | gmt pstext -R -J -DJ1c/1c -F+cTR -O -K -V >> $PS1
echo "Std = ${rmser_std} s" | gmt pstext -R -J -DJ-5.5c/2c -F+cTR+jLT -O -K -V >> $PS1
gmt pshistogram Deptd.txt -JX4.5i/3i -T0/40/2.5 -Bxf2.5a10+l"Depth (in km)" -By+l"Counts" -BWSne -W.5p -G0/128/128 -O -V -X5.8i -Y0.0i >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -E600 -Tf

evince $PS1


#gmt psbasemap -JX10i/6i -R0/200/0/60 -Bxf2.5a10+l"Depth (in km)" -Byf4a10+l"Vertical Error" -BWS+t"husn Depth Vs Vertical Error" -K -V >> Depxv.ps

#gmt psxy Depxv.txt -J -R -Sc0.2 -W0.5p,black -Gblue -O -V >> Depxv.ps

#gmt psbasemap -JX10i/6i -R0/200/0/60 -Bxf2.5a10+l"Depth (in km)" -Byf4a10+l"Horizontal Error" -BWS+t"husn Depth Vs Horiz. Error" -K -V >> Depxh.ps

#gmt psxy Depxh.txt -J -R -Sc0.2 -W0.5p,black -Gblue -O -V >> Depxh.ps

#gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

#gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tjef
#gmt psconvert Verr.ps -A+m3c/3c/3c/3c -Tg
#gmt psconvert Rmser.ps -A+m3c/3c/3c/3c -Tg
#gmt psconvert Deptd.ps -A+m3c/3c/3c/3c -Tg
#gmt psconvert Depxv.ps -A+m3c/3c/3c/3c -Tg
#gmt psconvert Depxh.ps -A+m3c/3c/3c/3c -Tg


#rm -f Depxy.txt Deptd.txt Rmser.txt Herr.txt Verr.txt test.txt test2.txt test3.txt Depxv.txt Depxh.txt
#EOF
