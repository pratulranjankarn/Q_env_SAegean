#!/bin/bash

rm -f ~/env_EGEL/intgrl_egel_compil_coda_nr.txt

for h in /home/user/env_EGEL/EGEL_events/*
do
cd $h
echo $h
#rm -f *.ascii
#rm -f out*.txt
rm -f env_integral_coda.txt

octave <<EOF

clc;
close all;

list = dir('out_N*envwin*.txt');

filnam = sprintf('env_integral_coda.txt');
fid = fopen(filnam,'w');

for i_iter = 1:length(list)
    
ay = strrep(list(i_iter).name,'out_N_','');

stnm = ay(1:4);
bnd = ay(end-6:end-4);
if (ay(4)=='_')
    stnm = ay(1:3);
end
if (bnd(1)=='-')
   bnd = ay(end-7:end-4);
   if (bnd(1)=='6')
	bnd = ay(end-8:end-4);
   end
end

temp = load(list(i_iter).name);

sdist = num2str(temp(end))

v_file=fullfile('/home/user/env_EGEL/station_files/',stnm,'vel_avg_husn.txt');
v_avg=load(v_file);
v_avg = v_avg/1000;

dta = temp(1:end-1);

if length(dta) > 10000
  fsamp = 100;
elseif (length(dta) > 5000) && (length(dta) < 10000)
  fsamp = 50;
  x = linspace(-1,100,length(dta));
  x_nw = linspace(-1,100,10100);
  dta_nw = interp1(x,dta,x_nw);
  dta = dta_nw;
else
  fsamp = 20;
  x = linspace(-1,100,length(dta));
  x_nw = linspace(-1,100,10100);
  dta_nw = interp1(x,dta,x_nw);
  dta = dta_nw;
end

fsamp=100;

samp_st = round(2*sdist/v_avg*fsamp);
samp_end = round(samp_st + 45*fsamp);

U4 = trapz(dta(fsamp*51:fsamp*61))/(10*fsamp);
U_C = trapz(dta(samp_st:samp_end))/U4;

U_C = num2str(U_C);

AB = {stnm, sdist, bnd, 'UC_norm', U_C};

fprintf(fid,"%s %s %s %s %s\n",AB'{:});

end
fclose(fid);
close all;
clear all;
quit
EOF
#SOD
cat env_integral_coda.txt >> /home/user/env_EGEL/intgrl_egel_compil_coda_nr.txt
cd ..
done
