clc;
close all;

list = dir('N_*_*HZ.sacii');

for i_iter = 1:length(list)

temp=replace(list(i_iter).name,'HZ','envwin');
temp1=replace(list(i_iter).name,'HZ','envwin_1-2');
temp2=replace(list(i_iter).name,'HZ','envwin_2-4');
temp3=replace(list(i_iter).name,'HZ','envwin_4-8');
temp4=replace(list(i_iter).name,'HZ','envwin_8-16');
temp5=replace(list(i_iter).name,'HZ','envwin_16-32');

temp1 = strcat('out_',replace(temp1,'sacii','txt'));
temp2 = strcat('out_',replace(temp2,'sacii','txt'));
temp3 = strcat('out_',replace(temp3,'sacii','txt'));
temp4 = strcat('out_',replace(temp4,'sacii','txt'));
temp5 = strcat('out_',replace(temp5,'sacii','txt'));

dta = load(temp1);

if length(dta) > 10000
fsamp = 100;
elseif (length(dta) > 5000) && (length(dta) < 10000)
fsamp = 50;
else
fsamp = 20;
end

wlen = round(2*fsamp+1);


figure(1);
subplot(5,1,1);
plot(dta(1:50*fsamp));
hold on;
plot(y(1:n_end));
title('1-2 Hz');

delete(temp1);

temp1 = replace(temp1,'.txt','_S.txt');
fid=fopen(temp1,'w');
fprintf(fid,'%f\n',y(1:parm));
fclose(fid);
temp1 = replace(temp1,'_S.txt','_Coda.txt');
fid=fopen(temp1,'w');
fprintf(fid,'%f\n',y(parm+1:n_end));
fclose(fid);

dta = load(temp2);
bw = hann(fsamp+1);
y = filtfilt(bw,1,dta);
dta = normc(dta);
y = normc(y);
r = snr(y,ones(size(y))*min(y));
r_set = 20 * log10((y.^2)./((min(y))^2));
n_f = find(r_set < (mean(r_set)/12));
parm = round(10.2 * fsamp);
n_end = min(n_f(n_f>parm))

subplot(5,1,2);
plot(dta(1:50*fsamp));
hold on;
plot(y(1:n_end));
title('2-4 Hz');

delete(temp2);

temp2 = replace(temp2,'.txt','_S.txt');
fid=fopen(temp2,'w');
fprintf(fid,'%f\n',y(1:parm));
fclose(fid);
temp2 = replace(temp2,'_S.txt','_Coda.txt');
fid=fopen(temp2,'w');
fprintf(fid,'%f\n',y(parm+1:n_end));
fclose(fid);

dta = load(temp3);
bw = hann(fsamp/2+1);
y = filtfilt(bw,1,dta);
dta = normc(dta);
y = normc(y);
r = snr(y,ones(size(y))*min(y));
r_set = 20 * log10((y.^2)./((min(y))^2));
n_f = find(r_set < (mean(r_set)/12));
parm = round(10.2 * fsamp);
n_end = min(n_f(n_f>parm))

subplot(5,1,3);
plot(dta(1:50*fsamp));
hold on;
plot(y(1:n_end));
title('4-8 Hz');

delete(temp3);

temp3 = replace(temp3,'.txt','_S.txt');
fid=fopen(temp3,'w');
fprintf(fid,'%f\n',y(1:parm));
fclose(fid);
temp3 = replace(temp3,'_S.txt','_Coda.txt');
fid=fopen(temp3,'w');
fprintf(fid,'%f\n',y(parm+1:n_end));
fclose(fid);

if fsamp > 20
dta = load(temp4);
bw = hann(round(fsamp/4+1));
y = filtfilt(bw,1,dta);
dta = normc(dta);
y = normc(y);
r = snr(y,ones(size(y))*min(y));
r_set = 20 * log10((y.^2)./((min(y))^2));
n_f = find(r_set < (mean(r_set)/12));
parm = round(10.2 * fsamp);
n_end = min(n_f(n_f>parm))

subplot(5,1,4);
plot(dta(1:50*fsamp));
hold on;
plot(y(1:n_end));
title('8-16 Hz');

delete(temp4);

temp4 = replace(temp4,'.txt','_S.txt');
fid=fopen(temp4,'w');
fprintf(fid,'%f\n',y(1:parm));
fclose(fid);
temp4 = replace(temp4,'_S.txt','_Coda.txt');
fid=fopen(temp4,'w');
fprintf(fid,'%f\n',y(parm+1:n_end));
fclose(fid);
end

if fsamp>50
dta = load(temp5);
bw = hann(round(fsamp/8+1));
y = filtfilt(bw,1,dta);
dta = normc(dta);
y = normc(y);
r = snr(y,ones(size(y))*min(y));
r_set = 20 * log10((y.^2)./((min(y))^2));
n_f = find(r_set < (mean(r_set)/12));
parm = round(10.2 * fsamp);
n_end = min(n_f(n_f>parm))

subplot(5,1,5);
plot(dta(1:50*fsamp));
hold on;
plot(y(1:n_end));
title('16-32 Hz');

delete(temp5);

temp5 = replace(temp5,'.txt','_S.txt');
fid=fopen(temp5,'w');
fprintf(fid,'%f\n',y(1:parm));
fclose(fid);
temp5 = replace(temp5,'_S.txt','_Coda.txt');
fid=fopen(temp5,'w');
fprintf(fid,'%f\n',y(parm+1:n_end));
fclose(fid);
end

temp = replace(temp,'sacii','png');
saveas(gcf,temp);

clf;
end