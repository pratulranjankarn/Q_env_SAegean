clc;
close all;

crcmrad = 62;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('aniso_st.txt','r');
% fid = fopen('Qsc12_final.txt','r');
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

ar_tri = zeros(1,length(lst));

lat = dt.Points(:,2);
lon = dt.Points(:,1);
h = zeros(length(lat),1);
ellp = referenceEllipsoid("wgs84");
[st_x,st_y,st_z] = geodetic2ecef(ellp,lat,lon,h);
for i=1:length(ar_tri)
tri_pt1 = [st_x(lst(i,1)) st_y(lst(i,1))];
tri_pt2 = [st_x(lst(i,2)) st_y(lst(i,2))];
tri_pt3 = [st_x(lst(i,3)) st_y(lst(i,3))];
ar_tri(i) = 0.5 * abs((tri_pt2(1)-tri_pt1(1))*(tri_pt3(2)-tri_pt1(2))-(tri_pt3(1)-tri_pt1(1))*(tri_pt2(2)-tri_pt1(2)));
end
% histogram(ar_tri)
idx = ar_tri>1.4*10^9;

% [crc,crr] = circumcenter(dt);
% crr = crr*111.1949;
% idx = find(crr > crcmrad);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));
triplot(dt2);

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

% fid = fopen('Qsc12_intrp.txt','w');
fid = fopen('aniso_st_intrp.txt','r');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qsc24_final.txt','r');
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
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

fid = fopen('Qsc24_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qsc48_final.txt','r');
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

fid = fopen('Qsc48_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qsc816_final.txt','r');
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

fid = fopen('Qsc816_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qi12_final.txt','r');
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

fid = fopen('Qi12_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qi24_final.txt','r');
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

fid = fopen('Qi24_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qi48_final.txt','r');
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

fid = fopen('Qi48_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Qi816_final.txt','r');
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

fid = fopen('Qi816_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);
