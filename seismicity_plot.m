close all;
clearvars;

wkdir = '/home/user/env_EGEL/env_gauss_7.2,3.6,1.8,.90_pkdbl';


abs_pdir = sprintf('%s/ray_prmtrs12.txt',wkdir);
fid_rp = fopen(abs_pdir,'r');
A = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/ray_prmtrs24.txt',wkdir);
fid_rp = fopen(abs_pdir,'r');
B = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/ray_prmtrs48.txt',wkdir);
fid_rp = fopen(abs_pdir,'r');
C = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/ray_prmtrs816.txt',wkdir);
fid_rp = fopen(abs_pdir,'r');
D = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);


mat1 = [A{2} A{1} A{3};
B{2} B{1} B{3};
C{2} C{1} C{3};
D{2} D{1} D{3}];


un_seis = unique(mat1,'rows','stable');


idx = find(un_seis(:,3)>100);

un_seis(idx,:) = [];

fid = fopen('seismicity_final.txt','w');
fprintf(fid,'%.6f %.6f %.6f\n',un_seis');
fclose(fid);

mat12 = [A{2} A{1} A{5} A{4}];

fid = fopen('ray_1-2.txt','w');

for i=1:length(mat12)
   fprintf(fid,'>-Z\n');
   fprintf(fid,'%f %f\n',mat12(i,3),mat12(i,4));
   fprintf(fid,'%f %f\n',mat12(i,1),mat12(i,2));
end

fclose(fid);

mat24 = [B{2} B{1} B{5} B{4}];

fid = fopen('ray_2-4.txt','w');

for i=1:length(mat24)
   fprintf(fid,'>-Z\n');
   fprintf(fid,'%f %f\n',mat24(i,3),mat24(i,4));
   fprintf(fid,'%f %f\n',mat24(i,1),mat24(i,2));
end

fclose(fid);

mat48 = [C{2} C{1} C{5} C{4}];

fid = fopen('ray_4-8.txt','w');

for i=1:length(mat48)
   fprintf(fid,'>-Z\n');
   fprintf(fid,'%f %f\n',mat48(i,3),mat48(i,4));
   fprintf(fid,'%f %f\n',mat48(i,1),mat48(i,2));
end

fclose(fid);

mat816 = [D{2} D{1} D{5} D{4}];

fid = fopen('ray_8-16.txt','w');

for i=1:length(mat816)
   fprintf(fid,'>-Z\n');
   fprintf(fid,'%f %f\n',mat816(i,3),mat816(i,4));
   fprintf(fid,'%f %f\n',mat816(i,1),mat816(i,2));
end

fclose(fid);