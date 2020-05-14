clc;
close all;

fid = fopen('Qc24_final.txt','r');
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
idx = find(crr > 50);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

Vq = Vq*1000;

ptset=[M(:,1) M(:,2) Vq];

fid = fopen('Qc24_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

fid = fopen('Qc48_final.txt','r');
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
idx = find(crr > 50);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

Vq = Vq*1000;

ptset=[M(:,1) M(:,2) Vq];

fid = fopen('Qc48_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

fid = fopen('Qc816_final.txt','r');
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
idx = find(crr > 50);
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = find(isnan(ti));
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

Vq = Vq*1000;

ptset=[M(:,1) M(:,2) Vq];

fid = fopen('Qc816_intrp.txt','w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);