#!/bin/bash

export HDF5_DISABLE_VERSION_CHECK=2

gmt set MAP_FRAME_PEN = thickest,black

gmt set FONT_LABEL = 14p,Helvetica,black

gmt set FONT_ANNOT_PRIMARY = 14p,Helvetica-Bold,black

gmt set MAP_FRAME_TYPE = fancy

gmt set PS_PAGE_ORIENTATION Landscape

PS1=geomap_SAegean2.ps
REGION="20/29.5/34.0/38.6"
SCALE=1.0i
GRID="/home/user/Downloads/ETOPO1_Bed_g_gmt4.grd/data"
INT="./topo.int"
CPT="custom.cpt"
seis_data="test.txt"

BLUE="65/105/225"

rm -f test.txt

#gmt makecpt -Cglobe -T-20000/20000/250 > $CPT
#gmt makecpt -Chaxby -T-5000/5000/25 > $CPT
rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt1.0 -nb+c+t1.0 -G$INT

gmt grdimage $GRID -I$INT -Ctopo_hillshade.cpt -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K > $PS1

gmt pscoast -Jm$SCALE -R$REGION -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/home/user/GSHHG/gshhg-gmt-2.3.7 -Dh -V -O -K >> $PS1

gmt psxy active_faults.txt -J -R -W1.2p,orange -V -O -K >> $PS1
#99/5/22

gmt psxy volcanoes.xy -J -R -W.7p,black -St0.5 -Gyellow -V -O -K >> $PS1
#207/16/32

#gmt psxy stlist_Qi816.txt -J -R -W.3p,black -Ss0.45 -Gyellow -V -O -K >> $PS1

#gmt psxy stlist_husn.xy -J -R -W.3p,black -Sc0.35 -Gred -V -O -K >> $PS1

gmt set FONT_ANNOT_PRIMARY = 15p,Helvetica,black

gmt set FONT_LABEL = 15p,Helvetica,black

gmt pslegend -DjLB+w2.41i/1.25i -F+gwhite -J -R -O -V << EOF >> $PS1
H 16,Helvetica-Bold Legend
S 0.1i t 0.50 yellow 0.7p,black 0.22i Active Volcanoes
G 0.03i
S 0.25i - 0.5i - 1.2p,orange 0.55i Active faults
G 0.03i
EOF

#gmt psconvert $PS1 -A+m3c/3c/3c/3c -E600 -Tg

gmt psconvert $PS1 -A+m3c/3c/3c/3c -E600 -Tf

evince $PS1


# S 0.1i c 0.35 red 0.3p,black 0.2i HUSN stations
# G 0.03i
# S 0.1i s 0.45 yellow 0.3p,black 0.2i EGELADOS stations
# G 0.03i