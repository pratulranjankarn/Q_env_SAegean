clc;
close all;
clearvars;

crcmrad = 62;

wkdir = '/home/user/env_EGEL/env_avg_chk2';

fid = fopen('stlist_f.txt','r');
st_db = textscan(fid,'%s %f %f');
fclose(fid);

ry_lmt = 10;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1-2 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_pdir = sprintf('%s/env_R_1-2.txt',wkdir);
fidR = fopen(abs_pdir,'r');
C = textscan(fidR,'%f %f %f');
fclose(fidR);

idcs_ev = unique(C{3});
ido_ttl = arrayfun(@(x)find(C{3}==x),idcs_ev,'UniformOutput',false);

Len = cellfun('length', ido_ttl);
[~,mxidx] = maxk(Len,3);

evmxst = mxidx(2);

set1 = C{1}(ido_ttl{evmxst});

idcs_st = unique(C{1});
idx = arrayfun(@(x)find(idcs_st==x),set1,'UniformOutput',false);

ido_ttl2 = arrayfun(@(x)find(C{1}==x),idcs_st,'UniformOutput',false);

set3 = cellfun(@(x)ido_ttl2{x},idx,'UniformOutput',false);



ste_efc = zeros(1,length(idcs));
for i=1:length(ido_ttl)
    s_g = length(ido_ttl{i});
    if (s_g > ry_lmt)
        tempg = C{2}(ido_ttl{i});
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        %%{
        for k = 1:100
            %k
            diffg = abs(tempg - mean_g);
            qw = find(diffg < 1.345*mad_g);
            sw = find(diffg >= 1.345*mad_g);
            tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
            mean_g = mean(tempg);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        ste_efc(i) = mean_g;       
    else
        ste_efc(i) = nan;
    end
end

aw = ~isnan(ste_efc);

abs_pdir = sprintf('%s/R_f_1-2.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',[st_db{3}(idcs(aw)) st_db{2}(idcs(aw)) ste_efc(aw)']');
fclose(fid);

fid = fopen(abs_pdir,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
[X,Y] = meshgrid(x,y);
c = cat(2,X',Y');
M_i = reshape(c,[],2);

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(M_i(:,1),M_i(:,2),msk_bx{1},msk_bx{2});

M = M_i(In_msk,:);

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[~,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = crr > crcmrad;
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = isnan(ti);
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

abs_pdir = sprintf('%s/R_f_1-2intrp.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2-4 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_pdir = sprintf('%s/env_R_2-4.txt',wkdir);
fidR = fopen(abs_pdir,'r');
C = textscan(fidR,'%f %f %f');
fclose(fidR);

idcs = unique(C{1});
ido_ttl = arrayfun(@(x)find(C{1}==x),idcs,'UniformOutput',false);

ste_efc = zeros(1,length(idcs));
for i=1:length(ido_ttl)
    s_g = length(ido_ttl{i});
    if (s_g > ry_lmt)
        tempg = C{2}(ido_ttl{i});
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        %%{
        for k = 1:100
            %k
            diffg = abs(tempg - mean_g);
            qw = find(diffg < 1.345*mad_g);
            sw = find(diffg >= 1.345*mad_g);
            tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
            mean_g = mean(tempg);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        ste_efc(i) = mean_g;       
    else
        ste_efc(i) = nan;
    end
end

aw = ~isnan(ste_efc);

abs_pdir = sprintf('%s/R_f_2-4.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',[st_db{3}(idcs(aw)) st_db{2}(idcs(aw)) ste_efc(aw)']');
fclose(fid);

fid = fopen(abs_pdir,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
[X,Y] = meshgrid(x,y);
c = cat(2,X',Y');
M_i = reshape(c,[],2);

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(M_i(:,1),M_i(:,2),msk_bx{1},msk_bx{2});

M = M_i(In_msk,:);

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
[~,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = crr > crcmrad;
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = isnan(ti);
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

abs_pdir = sprintf('%s/R_f_2-4intrp.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4-8 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_pdir = sprintf('%s/env_R_4-8.txt',wkdir);
fidR = fopen(abs_pdir,'r');
C = textscan(fidR,'%f %f %f');
fclose(fidR);

idcs = unique(C{1});
ido_ttl = arrayfun(@(x)find(C{1}==x),idcs,'UniformOutput',false);

ste_efc = zeros(1,length(idcs));
for i=1:length(ido_ttl)
    s_g = length(ido_ttl{i});
    if (s_g > ry_lmt)
        tempg = C{2}(ido_ttl{i});
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        %%{
        for k = 1:100
            %k
            diffg = abs(tempg - mean_g);
            qw = find(diffg < 1.345*mad_g);
            sw = find(diffg >= 1.345*mad_g);
            tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
            mean_g = mean(tempg);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        ste_efc(i) = mean_g;       
    else
        ste_efc(i) = nan;
    end
end

aw = ~isnan(ste_efc);

abs_pdir = sprintf('%s/R_f_4-8.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',[st_db{3}(idcs(aw)) st_db{2}(idcs(aw)) ste_efc(aw)']');
fclose(fid);


fid = fopen(abs_pdir,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
[X,Y] = meshgrid(x,y);
c = cat(2,X',Y');
M_i = reshape(c,[],2);

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(M_i(:,1),M_i(:,2),msk_bx{1},msk_bx{2});

M = M_i(In_msk,:);

dt = delaunayTriangulation(A{1},A{2});


lst = dt.ConnectivityList;
[~,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = crr > crcmrad;
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = isnan(ti);
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];

abs_pdir = sprintf('%s/R_f_4-8intrp.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8-16 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_pdir = sprintf('%s/env_R_8-16.txt',wkdir);
fidR = fopen(abs_pdir,'r');
C = textscan(fidR,'%f %f %f');
fclose(fidR);

idcs = unique(C{1});
ido_ttl = arrayfun(@(x)find(C{1}==x),idcs,'UniformOutput',false);

ste_efc = zeros(1,length(idcs));
for i=1:length(ido_ttl)
    s_g = length(ido_ttl{i});
    if (s_g > ry_lmt)
        tempg = C{2}(ido_ttl{i});
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        %%{
        for k = 1:100
            %k
            diffg = abs(tempg - mean_g);
            qw = find(diffg < 1.345*mad_g);
            sw = find(diffg >= 1.345*mad_g);
            tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
            mean_g = mean(tempg);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        ste_efc(i) = mean_g;       
    else
        ste_efc(i) = nan;
    end
end

aw = ~isnan(ste_efc);

abs_pdir = sprintf('%s/R_f_8-16.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',[st_db{3}(idcs(aw)) st_db{2}(idcs(aw)) ste_efc(aw)']');
fclose(fid);

fid = fopen(abs_pdir,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
[X,Y] = meshgrid(x,y);
c = cat(2,X',Y');
M_i = reshape(c,[],2);

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(M_i(:,1),M_i(:,2),msk_bx{1},msk_bx{2});

M = M_i(In_msk,:);

dt = delaunayTriangulation(A{1},A{2});

lst = dt.ConnectivityList;
[~,crr] = circumcenter(dt);
crr = crr*111.1949;
idx = crr > crcmrad;
lst(idx,:) = [];
dt2 = triangulation(lst,dt.Points);
%dt2 = delaunayTriangulation(pt(:,1),pt(:,2));

[ti,~] = pointLocation(dt2,M(:,1:2));
aw = isnan(ti);
M(aw,:) = [];
[ti,bc] = pointLocation(dt2,M(:,1:2));

triVals = A{3}(dt2(ti,:));
Vq = dot(bc',triVals')';

ptset=[M(:,1) M(:,2) Vq];


abs_pdir = sprintf('%s/R_f_8-16intrp.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);


%%
%%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16-32 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_pdir = sprintf('%s/env_R_16-32.txt',wkdir);
fidR = fopen(abs_pdir,'r');
C = textscan(fidR,'%f %f %f');
fclose(fidR);

idcs = unique(C{1});
ido_ttl = arrayfun(@(x)find(C{1}==x),idcs,'UniformOutput',false);

ste_efc = zeros(1,length(idcs));
for i=1:length(ido_ttl)
    s_g = length(ido_ttl{i});
    if (s_g > ry_lmt)
        tempg = C{2}(ido_ttl{i});
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        %%{
        for k = 1:100
            %k
            diffg = abs(tempg - mean_g);
            qw = find(diffg < 1.345*mad_g);
            sw = find(diffg >= 1.345*mad_g);
            tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
            mean_g = mean(tempg);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        ste_efc(i) = mean_g;       
    else
        ste_efc(i) = nan;
    end
end

aw = ~isnan(ste_efc);

abs_pdir = sprintf('%s/R_f_16-32.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',[st_db{3}(idcs(aw)) st_db{2}(idcs(aw)) ste_efc(aw)']');
fclose(fid);

fid = fopen(abs_pdir,'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

x=linspace(21,29.5,200);
y=linspace(34.5,38.5,200);
[X,Y] = meshgrid(x,y);
c = cat(2,X',Y');
M_i = reshape(c,[],2);

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(M_i(:,1),M_i(:,2),msk_bx{1},msk_bx{2});

M = M_i(In_msk,:);

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

abs_pdir = sprintf('%s/R_f_16-32intrp.txt',wkdir);
fid = fopen(abs_pdir,'w');
fprintf(fid,'%f %f %f\n',ptset');
fclose(fid);