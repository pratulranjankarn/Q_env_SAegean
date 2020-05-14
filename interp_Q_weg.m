clc;
close all;
clearvars;

wrk_path = '/home/user/env_EGEL/env_weg_test';

mkdir(wrk_path);

crcmrad = 62;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15g_1-2.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qsc12_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15g_2-4.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

%{

C = convexHull(dt);

factr = 50/111.1949/sqrt(2);

conv_hl = dt.Points(C,:);
pt = dt.Points;

ck_pt(1) = mean(conv_hl(:,1));
ck_pt(2) = mean(conv_hl(:,2));

for i=1:length(conv_hl)
if (conv_hl(i,1) < ck_pt(1)) && (conv_hl(i,2) < ck_pt(2))
conv_hl(i,:) = conv_hl(i,:) - 0.31796;
elseif (conv_hl(i,1) < ck_pt(1)) && (conv_hl(i,2) > ck_pt(2))
conv_hl(i,1) = conv_hl(i,1) - 0.31796;
conv_hl(i,2) = conv_hl(i,2) + 0.31796;
elseif (conv_hl(i,1) > ck_pt(1)) && (conv_hl(i,2) > ck_pt(2))
conv_hl(i,:) = conv_hl(i,:) + 0.31796;
else
conv_hl(i,1) = conv_hl(i,1) + 0.31796;
conv_hl(i,2) = conv_hl(i,2) - 0.31796;
end
end

pt(C,:) = conv_hl;

%}
lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
triplot(dt2);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qsc24_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15g_4-8.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});


lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qsc48_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15g_8-16.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qsc816_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15g_16-32.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qsc1632_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15b_1-2.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qi12_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15b_2-4.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qi24_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15b_4-8.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qi48_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15b_8-16.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qi816_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st1 = sprintf('%s/envnodes15b_16-32.txt',wrk_path);

fid = fopen(st1,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
M=[];
k=1;
for i=1:length(x)
for j=1:length(y)
M(k,:) = [x(i),y(j)];
k=k+1;
end
end

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[crc,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

st1 = sprintf('%s/Qi1632_intrp.txt',wrk_path);

fid = fopen(st1,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);
