
%%{
clc;
close all;
clearvars;

wrk_path = '/home/user/env_EGEL/env_avg_chk2';

mkdir(wrk_path);

%2-4 Hz 10 rays = 0.00012631
%2-4 Hz 15 rays = 0.00018947

%%%%%%%%% 1s intgrl 60km Hypo. distance 20 x 20 x 1
%{
sg_lim12 = 5;
sg_lim24 = 10;
sg_lim48 = 64;
sg_lim816 = 220;
sg_lim1632 = 159;
%}

%%%%%%%%% 1s intgrl 60km Hypo. distance 10 x 10 x 1
%{
sg_lim12 = 5;
sg_lim24 = 10;
sg_lim48 = 64;
sg_lim816 = 220;
sg_lim1632 = 159;
%}

%%%%%%%%% 2s intgrl 100km Hypo. distance 20 x 20 x 1
%{
sg_lim12 = 5;
sg_lim24 = 24;
sg_lim48 = 153;
sg_lim816 = 486;
sg_lim1632 = 320;
%}

%%%%%%%%% 2s intgrl 60km Hypo. distance 20 x 20 x 1
%{
sg_lim12 = 5;
sg_lim24 = 10;
sg_lim48 = 60;
sg_lim816 = 208;
sg_lim1632 = 162;
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 - 2 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
% fid = fopen('comb_12giwr_intgrl.txt','r');
% A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
% fclose(fid);

wgs84 = wgs84Ellipsoid('kilometers');

lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 40.0;
np_la = 31;
np_lo = 46;
np_de = 1;
Lat = linspace(latlow,lathig,np_la);
Long = linspace(lonlow,lonhig,np_lo);
%Depth = [45.0000 35.0000 25.0000 15.0000 5.0000];
Depth = [15.0000];
bincord_i = zeros(np_la*np_lo*np_de,3);
m=1;
for i = 1:np_la
for j =1: np_lo
for k=1:np_de
bincord_i(m,:) = [Lat(i),Long(j),Depth(k)];
m = m + 1;
end
end
end

bincord_i = single(bincord_i);

bincord = bincord_i;

%{

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

bincord_ii = bincord_i(In_msk,:);

%%{

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);

fid = fopen('NAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_iii(:,2),bincord_iii(:,1),msk_bx{1},msk_bx{2});

bincord = bincord_iii(~In_msk,:);

%}

n_bn = length(bincord);

[Nx, Ny, Nz] = geodetic2ecef(wgs84,bincord(:,1),bincord(:,2),bincord(:,3));

Nx = single(Nx);
Ny = single(Ny);
Nz = single(Nz);

options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
options.Display = 'off';

%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_12giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs12.txt',wrk_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

freq = 1.5;

P_bin12 = cell(length(A{1}),1);
g_ry = cell(length(A{1}),1);
b_ry = cell(length(A{1}),1);

sm_n= zeros(length(A{1}),1,'single');

factr = 10000000;

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) < 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2ecef(wgs84,A{1}(i),A{2}(i),A{3}(i));
   [Sx,Sy,Sz] = geodetic2ecef(wgs84,A{4}(i),A{5}(i),0);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
   dist = (Nx-Cx).^2 + (Ny-Cy).^2 + (Nz-Cz).^2;
   id = find(dist<(sm_n(i).^2));

   id15 = find(bincord(id,3)==15);
   
   temp_rw = zeros(1,n_bn,'uint8');
   
   temp_rw(id(id15)) = 1;
   
   P_bin12{i} = temp_rw;
   g_ry{i} = single(A{8}(i)*v_avg*factr);
   b_ry{i} = single(A{9}(i)*factr);
   else
      P_bin12{i} = uint8(P_bin12{i});
      g_ry{i} = single(g_ry{i});
      b_ry{i} = single(b_ry{i});
   end
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_1-2.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_1-2.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnode[In1,On] = inpolygon(bincord(:,2),bincord(:,1),msk_bx{1},msk_bx{2});s30b_1-2.txt','w');
%fid7 = fopen('envnodes35g_1-2.txt','w');
%fid8 = fopen('envnodes35b_1-2.txt','w');


P_bin12 = uint8(cell2mat(P_bin12));
g_ry = single(cell2mat(g_ry));
b_ry = single(cell2mat(b_ry));

% g_f = lsqlin(double(P_bin12),double((g_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);
% b_f = lsqlin(double(P_bin12),double((b_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);

g_f1 = single(P_bin12)\g_ry;
b_f1 = single(P_bin12)\b_ry;

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

g_f2 = g_f1(In_msk);
b_f2 = b_f1(In_msk);
bincord_ii = bincord_i(In_msk,:);

%%{

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);
g_f = g_f2(~In_msk);
b_f = b_f2(~In_msk);

%}

b_f = b_f/factr;
g_f = g_f/factr;

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

idx = find(Qsc_inv > 0 & Qi_inv >0);

un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qsc_inv(idx))];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qi_inv(idx))];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);


%%
%{
c_n= zeros(length(A{1}),2,'single');
r_n= zeros(length(A{1}),1,'single');
z_n= zeros(length(A{1}),1,'single');
sm_n= zeros(length(A{1}),1,'single');


%n_arr5 = cell(length(A{1}),1);
%n_arr10 = cell(length(A{1}),1);
n_arr15 = cell(length(A{1}),1);
%n_arr25 = cell(length(A{1}),1);
%n_arr30 = cell(length(A{1}),1);
%n_arr35 = cell(length(A{1}),1);

factr = 100000;

for i=1:length(A{1})
    %i
   v_avg = A{6}(i)/1000;
   waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   b = track(waypts);
   c_n(i,:) = b(15,:);
   r_n(i) = A{7}(i);
   sm_n(i) = (A{7}(i) + v_avg * A{11}(i))/2;
   snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Cx,Cy,Cz] = geodetic2ecef(wgs84,c_n(i,1),c_n(i,2),z_n(i));
   dist = (Nx-Cx).^2 + (Ny-Cy).^2 + (Nz-Cz).^2;
   id = find(dist<(sm_n(i).^2));
  
   %id5 = find(bincord(id,3)==5);
   %id10 = find(bincord(id,3)==10);
   id15 = find(bincord(id,3)==15);
   %id25 = find(bincord(id,3)==25);
   %id30 = find(bincord(id,3)==30);
   %id35 = find(bincord(id,3)==35);
   %n_arr5{i} = [uint32(id(id5)) uint32(ones(length(id5),1)*A{8}(i)*factr) uint32(ones(length(id5),1)*A{9}(i)*factr) uint32(ones(length(id5),1)*v_avg*factr)];
   %n_arr10{i} = [uint32(id(id10)) uint32(ones(length(id10),1)*A{8}(i)*factr) uint32(ones(length(id10),1)*A{9}(i)*factr) uint32(ones(length(id10),1)*v_avg*factr)];
   n_arr15{i} = [uint32(id(id15)) uint32(ones(length(id15),1)*A{8}(i)*factr) uint32(ones(length(id15),1)*A{9}(i)*factr) uint32(ones(length(id15),1)*v_avg*factr)];
   %n_arr25{i} = [uint32(id(id25)) uint32(ones(length(id25),1)*A{8}(i)*factr) uint32(ones(length(id25),1)*A{9}(i)*factr) uint32(ones(length(id25),1)*v_avg*factr)];
   %n_arr30{i} = [uint32(id(id30)) uint32(ones(length(id30),1)*A{8}(i)*factr) uint32(ones(length(id30),1)*A{9}(i)*factr) uint32(ones(length(id30),1)*v_avg*factr)];
   %n_arr35{i} = [uint32(id(id35)) uint32(ones(length(id35),1)*A{8}(i)*factr) uint32(ones(length(id35),1)*A{9}(i)*factr) uint32(ones(length(id35),1)*v_avg*factr)];   
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_1-2.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_1-2.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnodes30b_1-2.txt','w');
%fid7 = fopen('envnodes35g_1-2.txt','w');
%fid8 = fopen('envnodes35b_1-2.txt','w');


n_arr = uint32(cell2mat(n_arr15));
length(n_arr)
un_12 = uint32(unique(n_arr(:,1),'rows','stable'));
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
ido_ttl = arrayfun(@(x)find(n_arr(:,1)==x),un_12,'UniformOutput',false);

sz = zeros(1,length(ido_ttl));

for i = 1:length(ido_ttl)
sz(i) = length(ido_ttl{i});
end

figure(1);
histogram(sz,length(ido_ttl));

for j=1:length(un_12)
    %j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim12)
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        mad_b = mad(tempb,1)/0.67;
        mean_b = mean(tempb);
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
        
        for k = 1:100
            %k
            diffb = abs(tempb - mean_b);
            qw = find(diffb < 1.345*mad_b);
            sw = find(diffb >= 1.345*mad_b);
            tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
            mean_b = mean(tempb);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        g_f(j) = mean_g;
        b_f(j) = mean_b;
    else
        g_f(j) = NaN;
        b_f(j) = NaN;
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw

n_arr = single(cell2mat(n_arr30));
un_12 = uint32(unique(n_arr(:,1),'rows','stable'));
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
ido_ttl = arrayfun(@(x)find(n_arr(:,1)==x),un_12,'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim12)
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        mad_b = mad(tempb,1)/0.67;
        mean_b = mean(tempb);
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
        
        for k = 1:100
            %k
            diffb = abs(tempb - mean_b);
            qw = find(diffb < 1.345*mad_b);
            sw = find(diffb >= 1.345*mad_b);
            tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
            mean_b = mean(tempb);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        g_f(j) = mean_g;
        b_f(j) = mean_b;
    else
        g_f(j) = NaN;
        b_f(j) = NaN;
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fidt3,'%f %f %f\n',un_arr');
fclose(fidt3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fidt4,'%f %f %f\n',un_arr');
fclose(fidt4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw



n_arr = single(cell2mat(n_arr25));
un_12 = uint32(unique(n_arr(:,1),'rows','stable'));
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
ido_ttl = arrayfun(@(x)find(n_arr(:,1)==x),un_12,'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim12)
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        mad_b = mad(tempb,1)/0.67;
        mean_b = mean(tempb);
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
        
        for k = 1:100
            %k
            diffb = abs(tempb - mean_b);
            qw = find(diffb < 1.345*mad_b);
            sw = find(diffb >= 1.345*mad_b);
            tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
            mean_b = mean(tempb);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        g_f(j) = mean_g;
        b_f(j) = mean_b;
    else
        g_f(j) = NaN;
        b_f(j) = NaN;
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid6,'%f %f %f\n',un_arr');
fclose(fid6);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw

n_arr = single(cell2mat(n_arr35));
un_12 = uint32(unique(n_arr(:,1),'rows','stable'));
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
ido_ttl = arrayfun(@(x)find(n_arr(:,1)==x),un_12,'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim12)
        mad_g = mad(tempg,1)/0.67;
        mean_g = mean(tempg);
        mad_b = mad(tempb,1)/0.67;
        mean_b = mean(tempb);
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
        
        for k = 1:100
            %k
            diffb = abs(tempb - mean_b);
            qw = find(diffb < 1.345*mad_b);
            sw = find(diffb >= 1.345*mad_b);
            tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
            mean_b = mean(tempb);
            if (length(qw)==s_g)
                break;
            end
        end
        %%}
        g_f(j) = mean_g;
        b_f(j) = mean_b;
    else
        g_f(j) = NaN;
        b_f(j) = NaN;
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid7,'%f %f %f\n',un_arr');
fclose(fid7);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid8,'%f %f %f\n',un_arr');
fclose(fid8);

%%{

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 - 4 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_12giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs24.txt',wrk_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

freq = 1.5;

P_bin12 = cell(length(A{1}),1);
g_ry = cell(length(A{1}),1);
b_ry = cell(length(A{1}),1);

sm_n= zeros(length(A{1}),1,'single');

factr = 10000000;

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) < 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2ecef(wgs84,A{1}(i),A{2}(i),A{3}(i));
   [Sx,Sy,Sz] = geodetic2ecef(wgs84,A{4}(i),A{5}(i),0);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
   dist = (Nx-Cx).^2 + (Ny-Cy).^2 + (Nz-Cz).^2;
   id = find(dist<(sm_n(i).^2));

   id15 = find(bincord(id,3)==15);
   
   temp_rw = zeros(1,n_bn,'uint8');
   
   temp_rw(id(id15)) = 1;
   
   P_bin12{i} = temp_rw;
   g_ry{i} = single(A{8}(i)*v_avg*factr);
   b_ry{i} = single(A{9}(i)*factr);
   else
      P_bin12{i} = uint8(P_bin12{i});
      g_ry{i} = single(g_ry{i});
      b_ry{i} = single(b_ry{i});
   end
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_2-4.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_2-4.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnode[In1,On] = inpolygon(bincord(:,2),bincord(:,1),msk_bx{1},msk_bx{2});s30b_1-2.txt','w');
%fid7 = fopen('envnodes35g_1-2.txt','w');
%fid8 = fopen('envnodes35b_1-2.txt','w');


P_bin12 = uint8(cell2mat(P_bin12));
g_ry = single(cell2mat(g_ry));
b_ry = single(cell2mat(b_ry));

% g_f = lsqlin(double(P_bin12),double((g_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);
% b_f = lsqlin(double(P_bin12),double((b_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);

g_f1 = single(P_bin12)\g_ry;
b_f1 = single(P_bin12)\b_ry;

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

g_f2 = g_f1(In_msk);
b_f2 = b_f1(In_msk);
bincord_ii = bincord_i(In_msk,:);

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);
g_f = g_f2(~In_msk);
b_f = b_f2(~In_msk);

b_f = b_f/factr;
g_f = g_f/factr;

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

idx = find(Qsc_inv > 0 & Qi_inv >0);

un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qsc_inv(idx))];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qi_inv(idx))];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 - 8 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_12giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs48.txt',wrk_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

freq = 1.5;

P_bin12 = cell(length(A{1}),1);
g_ry = cell(length(A{1}),1);
b_ry = cell(length(A{1}),1);

sm_n= zeros(length(A{1}),1,'single');

factr = 10000000;

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) < 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2ecef(wgs84,A{1}(i),A{2}(i),A{3}(i));
   [Sx,Sy,Sz] = geodetic2ecef(wgs84,A{4}(i),A{5}(i),0);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
   dist = (Nx-Cx).^2 + (Ny-Cy).^2 + (Nz-Cz).^2;
   id = find(dist<(sm_n(i).^2));

   id15 = find(bincord(id,3)==15);
   
   temp_rw = zeros(1,n_bn,'uint8');
   
   temp_rw(id(id15)) = 1;
   
   P_bin12{i} = temp_rw;
   g_ry{i} = single(A{8}(i)*v_avg*factr);
   b_ry{i} = single(A{9}(i)*factr);
   else
      P_bin12{i} = uint8(P_bin12{i});
      g_ry{i} = single(g_ry{i});
      b_ry{i} = single(b_ry{i});
   end
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_4-8.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_4-8.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnode[In1,On] = inpolygon(bincord(:,2),bincord(:,1),msk_bx{1},msk_bx{2});s30b_1-2.txt','w');
%fid7 = fopen('envnodes35g_1-2.txt','w');
%fid8 = fopen('envnodes35b_1-2.txt','w');


P_bin12 = uint8(cell2mat(P_bin12));
g_ry = single(cell2mat(g_ry));
b_ry = single(cell2mat(b_ry));

% g_f = lsqlin(double(P_bin12),double((g_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);
% b_f = lsqlin(double(P_bin12),double((b_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);

g_f1 = single(P_bin12)\g_ry;
b_f1 = single(P_bin12)\b_ry;

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

g_f2 = g_f1(In_msk);
b_f2 = b_f1(In_msk);
bincord_ii = bincord_i(In_msk,:);

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);
g_f = g_f2(~In_msk);
b_f = b_f2(~In_msk);

b_f = b_f/factr;
g_f = g_f/factr;

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

idx = find(Qsc_inv > 0 & Qi_inv >0);

un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qsc_inv(idx))];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qi_inv(idx))];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8 - 16 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_12giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs816.txt',wrk_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

freq = 1.5;

P_bin12 = cell(length(A{1}),1);
g_ry = cell(length(A{1}),1);
b_ry = cell(length(A{1}),1);

sm_n= zeros(length(A{1}),1,'single');

factr = 10000000;

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) < 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2ecef(wgs84,A{1}(i),A{2}(i),A{3}(i));
   [Sx,Sy,Sz] = geodetic2ecef(wgs84,A{4}(i),A{5}(i),0);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
   dist = (Nx-Cx).^2 + (Ny-Cy).^2 + (Nz-Cz).^2;
   id = find(dist<(sm_n(i).^2));

   id15 = find(bincord(id,3)==15);
   
   temp_rw = zeros(1,n_bn,'uint8');
   
   temp_rw(id(id15)) = 1;
   
   P_bin12{i} = temp_rw;
   g_ry{i} = single(A{8}(i)*v_avg*factr);
   b_ry{i} = single(A{9}(i)*factr);
   else
      P_bin12{i} = uint8(P_bin12{i});
      g_ry{i} = single(g_ry{i});
      b_ry{i} = single(b_ry{i});
   end
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_8-16.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_8-16.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnode[In1,On] = inpolygon(bincord(:,2),bincord(:,1),msk_bx{1},msk_bx{2});s30b_1-2.txt','w');
%fid7 = fopen('envnodes35g_1-2.txt','w');
%fid8 = fopen('envnodes35b_1-2.txt','w');


P_bin12 = uint8(cell2mat(P_bin12));
g_ry = single(cell2mat(g_ry));
b_ry = single(cell2mat(b_ry));

% g_f = lsqlin(double(P_bin12),double((g_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);
% b_f = lsqlin(double(P_bin12),double((b_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);

g_f1 = single(P_bin12)\g_ry;
b_f1 = single(P_bin12)\b_ry;

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

g_f2 = g_f1(In_msk);
b_f2 = b_f1(In_msk);
bincord_ii = bincord_i(In_msk,:);

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);
g_f = g_f2(~In_msk);
b_f = b_f2(~In_msk);

b_f = b_f/factr;
g_f = g_f/factr;

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

idx = find(Qsc_inv > 0 & Qi_inv >0);

un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qsc_inv(idx))];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qi_inv(idx))];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16 - 32 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_12giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs1632.txt',wrk_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

freq = 1.5;

P_bin12 = cell(length(A{1}),1);
g_ry = cell(length(A{1}),1);
b_ry = cell(length(A{1}),1);

sm_n= zeros(length(A{1}),1,'single');

factr = 10000000;

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) < 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2ecef(wgs84,A{1}(i),A{2}(i),A{3}(i));
   [Sx,Sy,Sz] = geodetic2ecef(wgs84,A{4}(i),A{5}(i),0);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
   dist = (Nx-Cx).^2 + (Ny-Cy).^2 + (Nz-Cz).^2;
   id = find(dist<(sm_n(i).^2));

   id15 = find(bincord(id,3)==15);
   
   temp_rw = zeros(1,n_bn,'uint8');
   
   temp_rw(id(id15)) = 1;
   
   P_bin12{i} = temp_rw;
   g_ry{i} = single(A{8}(i)*v_avg*factr);
   b_ry{i} = single(A{9}(i)*factr);
   else
      P_bin12{i} = uint8(P_bin12{i});
      g_ry{i} = single(g_ry{i});
      b_ry{i} = single(b_ry{i});
   end
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_16-32.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_16-32.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnode[In1,On] = inpolygon(bincord(:,2),bincord(:,1),msk_bx{1},msk_bx{2});s30b_1-2.txt','w');
%fid7 = fopen('envnodes35g_1-2.txt','w');
%fid8 = fopen('envnodes35b_1-2.txt','w');


P_bin12 = uint8(cell2mat(P_bin12));
g_ry = single(cell2mat(g_ry));
b_ry = single(cell2mat(b_ry));

% g_f = lsqlin(double(P_bin12),double((g_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);
% b_f = lsqlin(double(P_bin12),double((b_ry)),[],[],[],[],...
%     ones(n_bn,1)*1e-6,ones(n_bn,1)*10,[],options);

g_f1 = single(P_bin12)\g_ry;
b_f1 = single(P_bin12)\b_ry;

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

g_f2 = g_f1(In_msk);
b_f2 = b_f1(In_msk);
bincord_ii = bincord_i(In_msk,:);

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);
g_f = g_f2(~In_msk);
b_f = b_f2(~In_msk);

b_f = b_f/factr;
g_f = g_f/factr;

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

idx = find(Qsc_inv > 0 & Qi_inv >0);

un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qsc_inv(idx))];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord_iii(idx,1) bincord_iii(idx,2) log10(Qi_inv(idx))];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);