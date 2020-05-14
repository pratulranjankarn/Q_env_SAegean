%% Assign values to bins inside the scattering ellipsoid/spheroid of each ray
%%{
clc;
close all;
clearvars;

wrk_path = '/home/user/env_EGEL/env_gauss_3_2_1.5_.75';

mkdir(wrk_path);

factr = 10000000;

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
lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 40.0;
Lat = latlow:0.20:lathig;
Long = lonlow:0.20:lonhig;
[X,Y] = meshgrid(Lat,Long);
temp = cat(2,X',Y');
bw = reshape(temp,length(Lat)*length(Long),2);
%Depth = [45.0000 35.0000 25.0000 15.0000 5.0000];
Depth = 15.0000;
bincord_i = [bw ones(length(bw),1)*Depth];

bincord_i = single(bincord_i);

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

bincord_ii = bincord_i(In_msk,:);

bincord = bincord_ii;

wgs84 = wgs84Ellipsoid('kilometers');

[Nx, Ny, Nz] = geodetic2enu(bincord(:,1),bincord(:,2),0,34,18,0,wgs84);

Nx = single(Nx);
Ny = single(Ny);
Nz = single(Nz);

n_grd = length(Nx);

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

sg_lim12 = 10;

sm_n= zeros(length(A{1}),1,'single');


%n_arr5 = cell(length(A{1}),1);
%n_arr10 = cell(length(A{1}),1);
n_arr15 = cell(length(A{1}),1);
%n_arr25 = cell(length(A{1}),1);
%n_arr30 = cell(length(A{1}),1);
%n_arr35 = cell(length(A{1}),1);



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
   [Ex,Ey,Ez] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
   [Sx,Sy,Sz] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
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
   n_arr15{i} = [uint32(id(id15)) uint32(ones(length(id15),1)*A{8}(i)*v_avg*factr) uint32(ones(length(id15),1)*A{9}(i)*factr) uint32(ones(length(id15),1)*v_avg*factr)];
   %n_arr25{i} = [uint32(id(id25)) uint32(ones(length(id25),1)*A{8}(i)*factr) uint32(ones(length(id25),1)*A{9}(i)*factr) uint32(ones(length(id25),1)*v_avg*factr)];
   %n_arr30{i} = [uint32(id(id30)) uint32(ones(length(id30),1)*A{8}(i)*factr) uint32(ones(length(id30),1)*A{9}(i)*factr) uint32(ones(length(id30),1)*v_avg*factr)];
   %n_arr35{i} = [uint32(id(id35)) uint32(ones(length(id35),1)*A{8}(i)*factr) uint32(ones(length(id35),1)*A{9}(i)*factr) uint32(ones(length(id35),1)*v_avg*factr)];   
   else
      n_arr15{i} = uint32(n_arr15{i});
   end
end

%fid1 = fopen('envnodes5g_1-2.txt','w');
%fid2 = fopen('envnodes5b_1-2.txt','w');
%fidt1 = fopen('envnodes10g_1-2.txt','w');
%fidt2 = fopen('envnodes10b_1-2.txt','w');
st1 = sprintf('%s/envnodes15g_1-2.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_1-2.txt',wrk_path);
st3 = sprintf('%s/envnodes15a_1-2.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');
%fid5 = fopen('envnodes25g_1-2.txt','w');
%fid6 = fopen('envnodes25b_1-2.txt','w');
%fidt3 = fopen('envnodes30g_1-2.txt','w');
%fidt4 = fopen('envnode[In1,On] = inpolygon(bincord(:,2),bincord(:,1),msk_bx{1},msk_bx{2});s30b_1-2.txt','w');
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

figure(2);
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

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(albdo)];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);


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

freq = 3;

sg_lim24 = 10;

sm_n= zeros(length(A{1}),1,'single');


%n_arr5 = cell(length(A{1}),1);
%n_arr10 = cell(length(A{1}),1);
n_arr15 = cell(length(A{1}),1);
%n_arr25 = cell(length(A{1}),1);
%n_arr30 = cell(length(A{1}),1);
%n_arr35 = cell(length(A{1}),1);



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
   [Ex,Ey,Ez] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
   [Sx,Sy,Sz] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
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
   n_arr15{i} = [uint32(id(id15)) uint32(ones(length(id15),1)*A{8}(i)*v_avg*factr) uint32(ones(length(id15),1)*A{9}(i)*factr) uint32(ones(length(id15),1)*v_avg*factr)];
   %n_arr25{i} = [uint32(id(id25)) uint32(ones(length(id25),1)*A{8}(i)*factr) uint32(ones(length(id25),1)*A{9}(i)*factr) uint32(ones(length(id25),1)*v_avg*factr)];
   %n_arr30{i} = [uint32(id(id30)) uint32(ones(length(id30),1)*A{8}(i)*factr) uint32(ones(length(id30),1)*A{9}(i)*factr) uint32(ones(length(id30),1)*v_avg*factr)];
   %n_arr35{i} = [uint32(id(id35)) uint32(ones(length(id35),1)*A{8}(i)*factr) uint32(ones(length(id35),1)*A{9}(i)*factr) uint32(ones(length(id35),1)*v_avg*factr)];   
   else
       n_arr15{i} = uint32(n_arr15{i});
   end
end

%fid1 = fopen('envnodes5g_2-4.txt','w');
%fid2 = fopen('envnodes5b_2-4.txt','w');
%fidt1 = fopen('envnodes10g_2-4.txt','w');
%fidt2 = fopen('envnodes10b_2-4.txt','w');
st1 = sprintf('%s/envnodes15g_2-4.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_2-4.txt',wrk_path);
st3 = sprintf('%s/envnodes15a_2-4.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');
%fid5 = fopen('envnodes25g_2-4.txt','w');
%fid6 = fopen('envnodes25b_2-4.txt','w');
%fidt3 = fopen('envnodes30g_2-4.txt','w');
%fidt4 = fopen('envnodes30b_2-4.txt','w');
%fid7 = fopen('envnodes35g_2-4.txt','w');
%fid8 = fopen('envnodes35b_2-4.txt','w');


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

figure(2);
histogram(sz,length(ido_ttl));

for j=1:length(un_12)
    %j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim24)
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

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(albdo)];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);

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

freq = 6;

sg_lim48 = 10;

sm_n= zeros(length(A{1}),1,'single');


%n_arr5 = cell(length(A{1}),1);
%n_arr10 = cell(length(A{1}),1);
n_arr15 = cell(length(A{1}),1);
%n_arr25 = cell(length(A{1}),1);
%n_arr30 = cell(length(A{1}),1);
%n_arr35 = cell(length(A{1}),1);



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
   [Ex,Ey,Ez] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
   [Sx,Sy,Sz] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
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
   n_arr15{i} = [uint32(id(id15)) uint32(ones(length(id15),1)*A{8}(i)*v_avg*factr) uint32(ones(length(id15),1)*A{9}(i)*factr) uint32(ones(length(id15),1)*v_avg*factr)];
   %n_arr25{i} = [uint32(id(id25)) uint32(ones(length(id25),1)*A{8}(i)*factr) uint32(ones(length(id25),1)*A{9}(i)*factr) uint32(ones(length(id25),1)*v_avg*factr)];
   %n_arr30{i} = [uint32(id(id30)) uint32(ones(length(id30),1)*A{8}(i)*factr) uint32(ones(length(id30),1)*A{9}(i)*factr) uint32(ones(length(id30),1)*v_avg*factr)];
   %n_arr35{i} = [uint32(id(id35)) uint32(ones(length(id35),1)*A{8}(i)*factr) uint32(ones(length(id35),1)*A{9}(i)*factr) uint32(ones(length(id35),1)*v_avg*factr)];   
   else
      n_arr15{i} = uint32(n_arr15{i});
   end
end

%fid1 = fopen('envnodes5g_4-8.txt','w');
%fid2 = fopen('envnodes5b_4-8.txt','w');
%fidt1 = fopen('envnodes10g_4-8.txt','w');
%fidt2 = fopen('envnodes10b_4-8.txt','w');
st1 = sprintf('%s/envnodes15g_4-8.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_4-8.txt',wrk_path);
st3 = sprintf('%s/envnodes15a_4-8.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');
%fid5 = fopen('envnodes25g_4-8.txt','w');
%fid6 = fopen('envnodes25b_4-8.txt','w');
%fidt3 = fopen('envnodes30g_4-8.txt','w');
%fidt4 = fopen('envnodes30b_4-8.txt','w');
%fid7 = fopen('envnodes35g_4-8.txt','w');
%fid8 = fopen('envnodes35b_4-8.txt','w');


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

figure(2);
histogram(sz,length(ido_ttl));

for j=1:length(un_12)
    %j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim48)
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

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(albdo)];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);

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

freq = 12;

sg_lim816 = 10;

sm_n= zeros(length(A{1}),1,'single');


%n_arr5 = cell(length(A{1}),1);
%n_arr10 = cell(length(A{1}),1);
n_arr15 = cell(length(A{1}),1);
%n_arr25 = cell(length(A{1}),1);
%n_arr30 = cell(length(A{1}),1);
%n_arr35 = cell(length(A{1}),1);



erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) <= 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
   [Sx,Sy,Sz] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
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
   n_arr15{i} = [uint32(id(id15)) uint32(ones(length(id15),1)*A{8}(i)*v_avg*factr) uint32(ones(length(id15),1)*A{9}(i)*factr) uint32(ones(length(id15),1)*v_avg*factr)];
   %n_arr25{i} = [uint32(id(id25)) uint32(ones(length(id25),1)*A{8}(i)*factr) uint32(ones(length(id25),1)*A{9}(i)*factr) uint32(ones(length(id25),1)*v_avg*factr)];
   %n_arr30{i} = [uint32(id(id30)) uint32(ones(length(id30),1)*A{8}(i)*factr) uint32(ones(length(id30),1)*A{9}(i)*factr) uint32(ones(length(id30),1)*v_avg*factr)];
   %n_arr35{i} = [uint32(id(id35)) uint32(ones(length(id35),1)*A{8}(i)*factr) uint32(ones(length(id35),1)*A{9}(i)*factr) uint32(ones(length(id35),1)*v_avg*factr)];   
   else
      n_arr15{i} = uint32(n_arr15{i});
   end
end

%fid1 = fopen('envnodes5g_8-16.txt','w');
%fid2 = fopen('envnodes5b_8-16.txt','w');
%fidt1 = fopen('envnodes10g_8-16.txt','w');
%fidt2 = fopen('envnodes10b_8-16.txt','w');
st1 = sprintf('%s/envnodes15g_8-16.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_8-16.txt',wrk_path);
st3 = sprintf('%s/envnodes15a_8-16.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');
%fid5 = fopen('envnodes25g_8-16.txt','w');
%fid6 = fopen('envnodes25b_8-16.txt','w');
%fidt3 = fopen('envnodes30g_8-16.txt','w');
%fidt4 = fopen('envnodes30b_8-16.txt','w');
%fid7 = fopen('envnodes35g_8-16.txt','w');
%fid8 = fopen('envnodes35b_8-16.txt','w');


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

figure(2);
histogram(sz,length(ido_ttl));

for j=1:length(un_12)
    %j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim816)
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

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(albdo)];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);

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

freq = 24;

sg_lim1632 = 10;

sm_n= zeros(length(A{1}),1,'single');


%n_arr5 = cell(length(A{1}),1);
%n_arr10 = cell(length(A{1}),1);
n_arr15 = cell(length(A{1}),1);
%n_arr25 = cell(length(A{1}),1);
%n_arr30 = cell(length(A{1}),1);
%n_arr35 = cell(length(A{1}),1);



erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:length(A{1})
    
   if (A{3}(i) <= 40)
    %i
   v_avg = A{6}(i);
   %waypts = single([A{1}(i),A{2}(i); A{4}(i),A{5}(i)]);
   %b = track(waypts);
   %c_n(i,:) = b(15,:);
   %r_n(i) = A{7}(i);
   sm_n(i) = (r_n(i) + v_avg * A{7}(i))/2;
   %snteta = (distance(A{1}(i),A{2}(i),c_n(i,1),c_n(i,2))*111.1949)/(r_n(i)/2);
   %z_n(i) = -(A{3}(i) - cos(asin(snteta))*r_n(i)/2);
   [Ex,Ey,Ez] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
   [Sx,Sy,Sz] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
   Cx = (Sx + Ex)/2;
   Cy = (Sy + Ey)/2;
   Cz = (Sz + Ez)/2;
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
   n_arr15{i} = [uint32(id(id15)) uint32(ones(length(id15),1)*A{8}(i)*v_avg*factr) uint32(ones(length(id15),1)*A{9}(i)*factr) uint32(ones(length(id15),1)*v_avg*factr)];
   %n_arr25{i} = [uint32(id(id25)) uint32(ones(length(id25),1)*A{8}(i)*factr) uint32(ones(length(id25),1)*A{9}(i)*factr) uint32(ones(length(id25),1)*v_avg*factr)];
   %n_arr30{i} = [uint32(id(id30)) uint32(ones(length(id30),1)*A{8}(i)*factr) uint32(ones(length(id30),1)*A{9}(i)*factr) uint32(ones(length(id30),1)*v_avg*factr)];
   %n_arr35{i} = [uint32(id(id35)) uint32(ones(length(id35),1)*A{8}(i)*factr) uint32(ones(length(id35),1)*A{9}(i)*factr) uint32(ones(length(id35),1)*v_avg*factr)];   
   else
      n_arr15{i} = uint32(n_arr15{i});
   end
end

%fid1 = fopen('envnodes5g_16-32.txt','w');
%fid2 = fopen('envnodes5b_16-32.txt','w');
%fidt1 = fopen('envnodes10g_16-32.txt','w');
%fidt2 = fopen('envnodes10b_16-32.txt','w');
st1 = sprintf('%s/envnodes15g_16-32.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_16-32.txt',wrk_path);
st3 = sprintf('%s/envnodes15a_16-32.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');
%fid5 = fopen('envnodes25g_16-32.txt','w');
%fid6 = fopen('envnodes25b_16-32.txt','w');
%fidt3 = fopen('envnodes30g_16-32.txt','w');
%fidt4 = fopen('envnodes30b_16-32.txt','w');
%fid7 = fopen('envnodes35g_16-32.txt','w');
%fid8 = fopen('envnodes35b_16-32.txt','w');


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

figure(2);
histogram(sz,length(ido_ttl));

for j=1:length(un_12)
    %j
    ido = uint32(ido_ttl{j});
    v_trk(j) = n_arr(ido(1),4)/factr;
    %length(ido)
    tempg = single(n_arr(ido,2))/factr;
    tempb = single(n_arr(ido,3))/factr;
    s_g = length(tempg);
    if (s_g > sg_lim1632)
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

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv);

un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qsc_inv)];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(Qi_inv)];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(un_12,1) bincord(un_12,2) log10(albdo)];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);

clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A
