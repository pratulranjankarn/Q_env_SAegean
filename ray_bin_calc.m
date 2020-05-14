clc;
close all;
clearvars;
%%
fid = fopen('ray_prmtrs12.txt','r');
A12 = textscan(fid,"%f %f %f %f %f %f %f %f %f");
fclose(fid);
%%
fid = fopen('ray_prmtrs24.txt','r');
A24 = textscan(fid,"%f %f %f %f %f %f %f %f %f");
fclose(fid);
%%
fid = fopen('ray_prmtrs48.txt','r');
A48 = textscan(fid,"%f %f %f %f %f %f %f %f %f");
fclose(fid);
%%
fid = fopen('ray_prmtrs816.txt','r');
A816 = textscan(fid,"%f %f %f %f %f %f %f %f %f");
fclose(fid);
%%
fid = fopen('ray_prmtrs1632.txt','r');
A1632 = textscan(fid,"%f %f %f %f %f %f %f %f %f");
fclose(fid);

%% The grid dimensions for the mapping

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
bincord = zeros(np_la*np_lo*np_de,3);
m=1;
for i = 1:np_la
for j =1: np_lo
for k=1:np_de
bincord(m,:) = [Lat(i),Long(j),Depth(k)];
m = m + 1;
end
end
end

bincord = single(bincord);

wgs84 = wgs84Ellipsoid('kilometers');

[Nx, Ny, Nz] = geodetic2ecef(wgs84,bincord(:,1),bincord(:,2),bincord(:,3));

Nx = single(Nx);
Ny = single(Ny);
Nz = single(Nz);

%%

wkdir = '/home/user/env_EGEL/env_gauss2D_2';

freq = 1.5;

ri =1;

g_bin12 = zeros(length(A12{1}),length(Nx),'single');
b_bin12 = zeros(length(A12{1}),length(Nx),'single');

sum_Pd = zeros(length(Nx),1);


for i=1:length(A12{1})
   ev_lat = A12{1}(i);
   ev_lon = A12{2}(i);
   ev_dep = A12{3}(i);
   
   st_lat = A12{4}(i);
   st_lon = A12{5}(i);
   
   v_avg = A12{6}(i);
   t_end = A12{7}(i);
   
   g_st = A12{8}(i);
   b = A12{9}(i);
   
   [ev_x,ev_y,ev_z] = geodetic2ecef(wgs84,ev_lat,ev_lon,ev_dep);
   [st_x,st_y,st_z] = geodetic2ecef(wgs84,st_lat,st_lon,0);
   Cx = (st_x + ev_x)/2;
   Cy = (st_y + ev_y)/2;
   Cz = (st_z + ev_z)/2;
   
   
   r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
   
   if (r < 500)
       
       s_mj = (v_avg * t_end + r)/2;
       s_mn = sqrt(s_mj^2 - (r/2)^2);
       
       slp_ry = (st_y-ev_y)/(st_x-ev_x);
       dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
       slp_rypp = -1/slp_ry;
       dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
       
       %%{
       P_sgma = [s_mj^2/9 s_mn^2/9];
       P_mu = [0 0];
       P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
       sum_Pd = sum_Pd + P_d;
       %}
       %{
       foc1_dst = sqrt((ev_x-Cx)^2 + (ev_y-Cy)^2);
       foc2_dst = sqrt((st_x-Cx)^2 + (st_y-Cy)^2);
       sum_dst = sqrt((dst_rypp-foc1_dst).^2 + (dst_ry-0).^2) + sqrt((dst_rypp-foc2_dst).^2 + (dst_ry-0).^2);
       min_dst = foc1_dst + foc2_dst;
       %P_d = normpdf(sum_dst,min_dst,s_mj/3);
       %}
       %%{
       dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
       
       idx = find(dst < 1);
       
       idx_l = length(idx);
       
       g_bin12(ri,idx) = P_d(idx) * g_st/sum(P_d);
       
       b_bin12(ri,idx) = P_d(idx) * b/sum(P_d);
       
       ri = ri + 1;
       %}
       %{
       pos_idx = find(sum_dst < 2*s_mj);
       sum_Pd(pos_idx) = P_d(pos_idx) + sum_Pd(pos_idx);
       g_st_set = (g_st * P_d);
       b_set = (b * P_d);
       count(pos_idx) = count(pos_idx) + 1;
       g_avg(pos_idx) = (g_avg(pos_idx) + g_st_set(pos_idx));
       b_avg(pos_idx) = (b_avg(pos_idx) + b_set(pos_idx));
       %}
   end
end

g_bin12((ri+1):end,:) = [];
b_bin12((ri+1):end,:) = [];

sg_lim12 = 5;

g_f = zeros(length(Nx),1);
b_f = zeros(length(Nx),1);

for j=1:length(Nx)
    %j
    tempg = g_bin12(:,j);
    tempb = b_bin12(:,j);
    tempg(tempg==0) = [];
    tempb(tempb==0) = [];
    s_g = length(tempg)
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

%g_f = g_f./sum_Pd;
%b_f = b_f./sum_Pd;


%{
%ray_bin12 = normc(ray_bin12);

b_fnl = ray_bin12\b_ry;
g_fnl = ray_bin12\g_st_ry;

b_fnl(b_fnl < 0) = 0;
b_fnl(b_fnl > 5) = 0;
g_fnl(g_fnl < 0) = 0;

%}

%%{
arset = [bincord log10(g_f.*(3.6)./(2*pi*freq)) log10(b_f./(2*pi*freq))];

sp = sprintf('%s/bin_31_46_1_1-2Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f\n',arset');
fclose(fid);
%}
