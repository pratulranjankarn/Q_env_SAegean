%% Assign values to bins inside the scattering ellipsoid/spheroid of...
%% each ray based on a gaussian distribution centered on ray midpoint
%%{
clc;
close all;
clearvars;

fl_path = '/home/user/env_EGEL/env_gauss_7.2,3.6,1.8,.90_pkdbl';

wrk_path = '/home/user/env_EGEL/env_gauss_7.2,3.6,1.8,.90_pkdbl';

mkdir(wrk_path);

factr = 10000000;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

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
Lat = latlow:0.10:lathig;
Long = lonlow:0.10:lonhig;
[X,Y] = meshgrid(Lat,Long);
temp = cat(2,X',Y');
bw = reshape(temp,length(Lat)*length(Long),2);
%Depth = [45.0000 35.0000 25.0000 15.0000 5.0000];
Depth = 15.0000;
bincord_i = [bw ones(length(bw),1)*Depth];

bincord_i = single(bincord_i);

%%{

fid = fopen('SAegean_poly_coord.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_i(:,2),bincord_i(:,1),msk_bx{1},msk_bx{2});

bincord_ii = bincord_i(In_msk,:);

bincord = bincord_ii;

%{

fid = fopen('SAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord_iii = bincord_ii(~In_msk,:);
%}
%{
fid = fopen('NAegean_sea_mask.txt','r');
msk_bx = textscan(fid,'%f %f');
fclose(fid);

[In_msk,~] = inpolygon(bincord_ii(:,2),bincord_ii(:,1),msk_bx{1},msk_bx{2});

bincord = bincord_ii(~In_msk,:);

%}

[Nx, Ny, Nz] = geodetic2enu(bincord(:,1),bincord(:,2),bincord(:,3),34.0,18,0,wgs84);

Nx = single(Nx);
Ny = single(Ny);
Nz = single(Nz);

n_grd = length(Nx);

fctr_mjr = 1/9; 
fctr_mnr = 1/9;

sgmy='3_3';

%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_12giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs12.txt',fl_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

n_a = length(A{1});

freq = 1.5;

P_bin12 = zeros(n_a,n_grd,'single');

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:n_a 
   if (A{3}(i) < 40)% && (r_n(i) < 60)
       %i
       t_end = A{7}(i);
       r = r_n(i);
       v_avg = A{6}(i);
       [ev_x,ev_y,ev_z] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
       [st_x,st_y,st_z] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
       s_mj = (v_avg * t_end + r)/2;
       s_mn = sqrt(4*s_mj^2 - r^2)/2;
       Cx = (st_x + ev_x)/2;
       Cy = (st_y + ev_y)/2;
       %%{
       slp_ry = (st_y-ev_y)/(st_x-ev_x);
       dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
       slp_rypp = -1/slp_ry;
       dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
       
       P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
       P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
       %{
       P_d = normpdf(dst_ry,0,s_mj/3);
       
       dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
       
       idx = find(dst < 1);
       
       idx2 = find(dst_rypp > r/2 & dst < 1);
       
       P_bin12(i,idx) = single(P_d(idx));
       P_bin12(i,idx2) = single(normpdf(dst_rypp(idx2),r/2,s_mj/3));
       %}
       P_bin12(i,:) = single(P_d/sum(P_d));
   end
end

% figure(1); plot3([A{2}(i),A{5}(i)],[A{1}(i),A{4}(i)],[0,0]); hold on; scatter3(bincord(:,2),...
% bincord(:,1),dst_ry); hold off;
% figure(2); plot3([A{2}(i),A{5}(i)],[A{1}(i),A{4}(i)],[0,0]); hold on; scatter3(bincord(:,2),...
% bincord(:,1),dst_rypp); hold off;
% figure(4); plot3([A{2}(i),A{5}(i)],[A{1}(i),A{4}(i)],[0,0]); hold on; scatter3(bincord(:,2),...
% bincord(:,1),P_bin12(i,:)); hold off;

% f = waitbar(0,'1','Name','Iterating through parameters for 1-2 Hz',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

%{
gs_set = logspace(-7,1,100);
bs_set = logspace(-7,1,100);
err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,3);
gs_tmp = gs_set(idx);

gs_set = logspace(log10(min(gs_tmp)),log10(max(gs_tmp)),100);

[~,idx] = mink(err_bset,3);

bs_tmp = bs_set(idx);

bs_set = logspace(log10(min(bs_tmp)),log10(max(bs_tmp)),100);

err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,1);

g_fi = gs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*gs_set(idx);

[~,idx] = mink(err_bset,1);

b_fi = bs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*bs_set(idx);

%%{

std_g = single(0.0025);
std_b = single(0.0025);

g_fmt = P_bin12 * g_fi;
b_fmt = P_bin12 * b_fi;

g_f = g_fi;
b_f = b_fi;

L_b = (sum(del2(b_f).^2)/n_grd)*n_a;
L_g = (sum(del2(g_f).^2)/n_grd)*n_a;

err_g = (g_fmt-A{8})'*(g_fmt-A{8}) + L_g;
err_b = (b_fmt-A{9})'*(b_fmt-A{9}) + L_b;

g_ff = single(0);
b_ff = single(0);

for iter=1:250000
   var = single(randi(n_grd));
   g_f2 = g_f;
   b_f2 = b_f;
   chus = randi(2);
   if chus == 1
       temp1 = g_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_g,(1-temp1)/std_g));
       g_f2(var) = temp1 + temp2 * std_g;
       L_g2 = (sum(del2(g_f2).^2)/n_grd)*n_a;
       g_fmt2 = P_bin12 * g_f2;
       err_g2 = (g_fmt2-A{8})'*(g_fmt2-A{8}) + L_g2;
       accp_prob = min(1,exp(-(err_g2-err_g)));
       if (accp_prob == 1 || accp_prob > rand(1))
           g_f = g_f2;
           err_g = err_g2;
       end
   else
       temp1 = b_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_b,(1-temp1)/std_b));
       b_f2(var) = temp1 + temp2 * std_b;
       L_b2 = (sum(del2(b_f2).^2)/n_grd)*n_a;
       b_fmt2 = P_bin12 * b_f2;
       err_b2 = (b_fmt2-A{9})'*(b_fmt2-A{9}) + L_b2;
       accp_prob = min(1,exp(-(err_b2-err_b)));
       if (accp_prob == 1 || accp_prob > rand(1))
           b_f = b_f2;
           err_b = err_b2;
       end
   end
   if (iter > 200000)
      g_ff = (g_f + g_ff*(iter-200000-1))/(iter-200000);
      b_ff = (b_f + b_ff*(iter-200000-1))/(iter-200000);
      iter
   end
%    waitbar(iter/250000,f,sprintf('%f:%f:%f',iter,err_g,err_b));
%    if getappdata(f,'canceling')
%         break
%    end
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

g_f = g_ff;
b_f = b_ff;
%}

st1 = sprintf('%s/envnodes15g_1-2_%s.txt',wrk_path,sgmy);
st2 = sprintf('%s/envnodes15b_1-2_%s.txt',wrk_path,sgmy);
st3 = sprintf('%s/envnodes15a_1-2_%s.txt',wrk_path,sgmy);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv)*100;

un_arr = [bincord(:,1) bincord(:,2) log10(Qsc_inv')];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(:,1) bincord(:,2) log10(Qi_inv')];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(:,1) bincord(:,2) albdo'];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);


%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

%fid = fopen('comb_24giwr_intgrl.txt','r');
%A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
%fclose(fid);

fp1 = sprintf('%s/ray_prmtrs24.txt',fl_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

n_a = length(A{1});

freq = 3;

P_bin12 = zeros(n_a,n_grd,'single');

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:n_a 
   if (A{3}(i) < 40)% && (r_n(i) < 60)
       %i
       t_end = A{7}(i);
       r = r_n(i);
       v_avg = A{6}(i);
       [ev_x,ev_y,ev_z] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
       [st_x,st_y,st_z] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
       s_mj = (v_avg * t_end + r)/2;
       s_mn = sqrt(4*s_mj^2 - r^2)/2;
       Cx = (st_x + ev_x)/2;
       Cy = (st_y + ev_y)/2;
       %%{
       slp_ry = (st_y-ev_y)/(st_x-ev_x);
       dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
       slp_rypp = -1/slp_ry;
       dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
       
       P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
       P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
       %{
       P_d = normpdf(dst_ry,0,s_mj/3);
       
       dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
       
       idx = find(dst < 1);
       
       idx2 = find(dst_rypp > r/2 & dst < 1);
       
       P_bin12(i,idx) = single(P_d(idx));
       P_bin12(i,idx2) = single(normpdf(dst_rypp(idx2),r/2,s_mj/3));
       %}
       P_bin12(i,:) = single(P_d/sum(P_d));
   end
end

g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

% f = waitbar(0,'1','Name','Iterating through parameters for 2-4 Hz',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%{
gs_set = logspace(-7,1,100);
bs_set = logspace(-7,1,100);
err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,3);
gs_tmp = gs_set(idx);

gs_set = logspace(log10(min(gs_tmp)),log10(max(gs_tmp)),100);

[~,idx] = mink(err_bset,3);

bs_tmp = bs_set(idx);

bs_set = logspace(log10(min(bs_tmp)),log10(max(bs_tmp)),100);

err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,1);

g_fi = single(gs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*gs_set(idx));

[~,idx] = mink(err_bset,1);

b_fi = single(bs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*bs_set(idx));

% g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
% b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

%%{

std_g = single(0.0025);
std_b = single(0.0025);

g_fmt = P_bin12 * g_fi;
b_fmt = P_bin12 * b_fi;

g_f = g_fi;
b_f = b_fi;

L_b = (sum(del2(b_f).^2)/n_grd)*n_a;
L_g = (sum(del2(g_f).^2)/n_grd)*n_a;

err_g = (g_fmt-A{8})'*(g_fmt-A{8}) + L_g;
err_b = (b_fmt-A{9})'*(b_fmt-A{9}) + L_b;

g_ff = single(0);
b_ff = single(0);

for iter=1:250000
   var = single(randi(n_grd));
   g_f2 = g_f;
   b_f2 = b_f;
   chus = randi(2);
   if chus == 1
       temp1 = g_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_g,(1-temp1)/std_g));
       g_f2(var) = temp1 + temp2 * std_g;
       L_g2 = (sum(del2(g_f2).^2)/n_grd)*n_a;
       g_fmt2 = P_bin12 * g_f2;
       err_g2 = (g_fmt2-A{8})'*(g_fmt2-A{8}) + L_g2;
       accp_prob = min(1,exp(-(err_g2-err_g)));
       if (accp_prob == 1 || accp_prob > rand(1))
           g_f = g_f2;
           err_g = err_g2;
       end
   else
       temp1 = b_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_b,(1-temp1)/std_b));
       b_f2(var) = temp1 + temp2 * std_b;
       L_b2 = (sum(del2(b_f2).^2)/n_grd)*n_a;
       b_fmt2 = P_bin12 * b_f2;
       err_b2 = (b_fmt2-A{9})'*(b_fmt2-A{9}) + L_b2;
       accp_prob = min(1,exp(-(err_b2-err_b)));
       if (accp_prob == 1 || accp_prob > rand(1))
           b_f = b_f2;
           err_b = err_b2;
       end
   end
   if (iter > 200000)
      g_ff = (g_f + g_ff*(iter-200000-1))/(iter-200000);
      b_ff = (b_f + b_ff*(iter-200000-1))/(iter-200000);
      iter
   end
%    waitbar(iter/250000,f,sprintf('%f:%f:%f',iter,err_g,err_b));
%    if getappdata(f,'canceling')
%         break
%    end
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

g_f = g_ff;
b_f = b_ff;
%}

st1 = sprintf('%s/envnodes15g_2-4_%s.txt',wrk_path,sgmy);
st2 = sprintf('%s/envnodes15b_2-4_%s.txt',wrk_path,sgmy);
st3 = sprintf('%s/envnodes15a_2-4_%s.txt',wrk_path,sgmy);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv)*100;

un_arr = [bincord(:,1) bincord(:,2) log10(Qsc_inv')];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(:,1) bincord(:,2) log10(Qi_inv')];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(:,1) bincord(:,2) albdo'];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 - 8 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

fp1 = sprintf('%s/ray_prmtrs48.txt',fl_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

n_a = length(A{1});

freq = 6;

P_bin12 = zeros(n_a,n_grd,'single');

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:n_a 
   if (A{3}(i) < 40)% && (r_n(i) < 60)
       %i
       t_end = A{7}(i);
       r = r_n(i);
       v_avg = A{6}(i);
       [ev_x,ev_y,ev_z] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
       [st_x,st_y,st_z] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
       s_mj = (v_avg * t_end + r)/2;
       s_mn = sqrt(4*s_mj^2 - r^2)/2;
       Cx = (st_x + ev_x)/2;
       Cy = (st_y + ev_y)/2;
       %%{
       slp_ry = (st_y-ev_y)/(st_x-ev_x);
       dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
       slp_rypp = -1/slp_ry;
       dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
       
       P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
       P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
       %{
       P_d = normpdf(dst_ry,0,s_mj/3);
       
       dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
       
       idx = find(dst < 1);
       
       idx2 = find(dst_rypp > r/2 & dst < 1);
       
       P_bin12(i,idx) = single(P_d(idx));
       P_bin12(i,idx2) = single(normpdf(dst_rypp(idx2),r/2,s_mj/3));
       %}
       P_bin12(i,:) = single(P_d/sum(P_d));
   end
end

g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

% f = waitbar(0,'1','Name','Iterating through parameters for 4-8 Hz',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%{
gs_set = logspace(-7,1,100);
bs_set = logspace(-7,1,100);
err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,3);
gs_tmp = gs_set(idx);

gs_set = logspace(log10(min(gs_tmp)),log10(max(gs_tmp)),100);

[~,idx] = mink(err_bset,3);

bs_tmp = bs_set(idx);

bs_set = logspace(log10(min(bs_tmp)),log10(max(bs_tmp)),100);

err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,1);

g_fi = single(gs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*gs_set(idx));

[~,idx] = mink(err_bset,1);

b_fi = single(bs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*bs_set(idx));

% g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
% b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

%%{

std_g = single(0.0025);
std_b = single(0.0025);

g_fmt = P_bin12 * g_fi;
b_fmt = P_bin12 * b_fi;

g_f = g_fi;
b_f = b_fi;

L_b = (sum(del2(b_f).^2)/n_grd)*n_a;
L_g = (sum(del2(g_f).^2)/n_grd)*n_a;

err_g = (g_fmt-A{8})'*(g_fmt-A{8}) + L_g;
err_b = (b_fmt-A{9})'*(b_fmt-A{9}) + L_b;

g_ff = single(0);
b_ff = single(0);

for iter=1:250000
   var = single(randi(n_grd));
   g_f2 = g_f;
   b_f2 = b_f;
   chus = randi(2);
   if chus == 1
       temp1 = g_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_g,(1-temp1)/std_g));
       g_f2(var) = temp1 + temp2 * std_g;
       L_g2 = (sum(del2(g_f2).^2)/n_grd)*n_a;
       g_fmt2 = P_bin12 * g_f2;
       err_g2 = (g_fmt2-A{8})'*(g_fmt2-A{8}) + L_g2;
       accp_prob = min(1,exp(-(err_g2-err_g)));
       if (accp_prob == 1 || accp_prob > rand(1))
           g_f = g_f2;
           err_g = err_g2;
       end
   else
       temp1 = b_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_b,(1-temp1)/std_b));
       b_f2(var) = temp1 + temp2 * std_b;
       L_b2 = (sum(del2(b_f2).^2)/n_grd)*n_a;
       b_fmt2 = P_bin12 * b_f2;
       err_b2 = (b_fmt2-A{9})'*(b_fmt2-A{9}) + L_b2;
       accp_prob = min(1,exp(-(err_b2-err_b)));
       if (accp_prob == 1 || accp_prob > rand(1))
           b_f = b_f2;
           err_b = err_b2;
       end
   end
   if (iter > 200000)
      g_ff = (g_f + g_ff*(iter-200000-1))/(iter-200000);
      b_ff = (b_f + b_ff*(iter-200000-1))/(iter-200000);
      iter
   end
%    waitbar(iter/250000,f,sprintf('%f:%f:%f',iter,err_g,err_b));
%    if getappdata(f,'canceling')
%         break
%    end
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

g_f = g_ff;
b_f = b_ff;
%}

st1 = sprintf('%s/envnodes15g_4-8_%s.txt',wrk_path,sgmy);
st2 = sprintf('%s/envnodes15b_4-8_%s.txt',wrk_path,sgmy);
st3 = sprintf('%s/envnodes15a_4-8_%s.txt',wrk_path,sgmy);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv)*100;

un_arr = [bincord(:,1) bincord(:,2) log10(Qsc_inv')];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(:,1) bincord(:,2) log10(Qi_inv')];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(:,1) bincord(:,2) albdo'];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8 - 16 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

fp1 = sprintf('%s/ray_prmtrs816.txt',fl_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

n_a = length(A{1});

freq = 12;

P_bin12 = zeros(n_a,n_grd,'single');

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:n_a 
   if (A{3}(i) < 40)% && (r_n(i) < 60)
       %i
       t_end = A{7}(i);
       r = r_n(i);
       v_avg = A{6}(i);
       [ev_x,ev_y,ev_z] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
       [st_x,st_y,st_z] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
       s_mj = (v_avg * t_end + r)/2;
       s_mn = sqrt(4*s_mj^2 - r^2)/2;
       Cx = (st_x + ev_x)/2;
       Cy = (st_y + ev_y)/2;
       %%{
       slp_ry = (st_y-ev_y)/(st_x-ev_x);
       dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
       slp_rypp = -1/slp_ry;
       dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
       
       P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
       P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
       %{
       P_d = normpdf(dst_ry,0,s_mj/3);
       
       dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
       
       idx = find(dst < 1);
       
       idx2 = find(dst_rypp > r/2 & dst < 1);
       
       P_bin12(i,idx) = single(P_d(idx));
       P_bin12(i,idx2) = single(normpdf(dst_rypp(idx2),r/2,s_mj/3));
       %}
       P_bin12(i,:) = single(P_d/sum(P_d));
   end
end

g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

% f = waitbar(0,'1','Name','Iterating through parameters for 8-16 Hz',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%{
gs_set = logspace(-7,1,100);
bs_set = logspace(-7,1,100);
err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,3);
gs_tmp = gs_set(idx);

gs_set = logspace(log10(min(gs_tmp)),log10(max(gs_tmp)),100);

[~,idx] = mink(err_bset,3);

bs_tmp = bs_set(idx);

bs_set = logspace(log10(min(bs_tmp)),log10(max(bs_tmp)),100);

err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,1);

g_fi = single(gs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*gs_set(idx));

[~,idx] = mink(err_bset,1);

b_fi = single(bs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*bs_set(idx));

% g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
% b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

%%{

std_g = single(0.0025);
std_b = single(0.0025);

g_fmt = P_bin12 * g_fi;
b_fmt = P_bin12 * b_fi;

g_f = g_fi;
b_f = b_fi;

L_b = (sum(del2(b_f).^2)/n_grd)*n_a;
L_g = (sum(del2(g_f).^2)/n_grd)*n_a;

err_g = (g_fmt-A{8})'*(g_fmt-A{8}) + L_g;
err_b = (b_fmt-A{9})'*(b_fmt-A{9}) + L_b;

g_ff = single(0);
b_ff = single(0);

for iter=1:250000
   var = single(randi(n_grd));
   g_f2 = g_f;
   b_f2 = b_f;
   chus = randi(2);
   if chus == 1
       temp1 = g_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_g,(1-temp1)/std_g));
       g_f2(var) = temp1 + temp2 * std_g;
       L_g2 = (sum(del2(g_f2).^2)/n_grd)*n_a;
       g_fmt2 = P_bin12 * g_f2;
       err_g2 = (g_fmt2-A{8})'*(g_fmt2-A{8}) + L_g2;
       accp_prob = min(1,exp(-(err_g2-err_g)));
       if (accp_prob == 1 || accp_prob > rand(1))
           g_f = g_f2;
           err_g = err_g2;
       end
   else
       temp1 = b_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_b,(1-temp1)/std_b));
       b_f2(var) = temp1 + temp2 * std_b;
       L_b2 = (sum(del2(b_f2).^2)/n_grd)*n_a;
       b_fmt2 = P_bin12 * b_f2;
       err_b2 = (b_fmt2-A{9})'*(b_fmt2-A{9}) + L_b2;
       accp_prob = min(1,exp(-(err_b2-err_b)));
       if (accp_prob == 1 || accp_prob > rand(1))
           b_f = b_f2;
           err_b = err_b2;
       end
   end
   if (iter > 200000)
      g_ff = (g_f + g_ff*(iter-200000-1))/(iter-200000);
      b_ff = (b_f + b_ff*(iter-200000-1))/(iter-200000);
      iter
   end
%    waitbar(iter/250000,f,sprintf('%f:%f:%f',iter,err_g,err_b));
%    if getappdata(f,'canceling')
%         break
%    end
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

g_f = g_ff;
b_f = b_ff;
%}

st1 = sprintf('%s/envnodes15g_8-16_%s.txt',wrk_path,sgmy);
st2 = sprintf('%s/envnodes15b_8-16_%s.txt',wrk_path,sgmy);
st3 = sprintf('%s/envnodes15a_8-16_%s.txt',wrk_path,sgmy);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv)*100;

un_arr = [bincord(:,1) bincord(:,2) log10(Qsc_inv')];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(:,1) bincord(:,2) log10(Qi_inv')];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(:,1) bincord(:,2) albdo'];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16 - 32 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear c_n r_n z_n sm_n n_arr5 n_arr15 n_arr25 n_arr35 A

fp1 = sprintf('%s/ray_prmtrs1632.txt',fl_path);
fid = fopen(fp1,'r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f');
fclose(fid);

n_a = length(A{1});

freq = 24;

P_bin12 = zeros(n_a,n_grd,'single');

erth_rad = 6371;

arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);

r_n = sqrt(arclen.^2 + A{3}.^2);

for i=1:n_a 
   if (A{3}(i) < 40)% && (r_n(i) < 60)
       %i
       t_end = A{7}(i);
       r = r_n(i);
       v_avg = A{6}(i);
       [ev_x,ev_y,ev_z] = geodetic2enu(A{1}(i),A{2}(i),A{3}(i),34.0,18,0,wgs84);
       [st_x,st_y,st_z] = geodetic2enu(A{4}(i),A{5}(i),0,34.0,18,0,wgs84);
       s_mj = (v_avg * t_end + r)/2;
       s_mn = sqrt(4*s_mj^2 - r^2)/2;
       Cx = (st_x + ev_x)/2;
       Cy = (st_y + ev_y)/2;
       %%{
       slp_ry = (st_y-ev_y)/(st_x-ev_x);
       dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
       slp_rypp = -1/slp_ry;
       dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
       
       P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
       P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
       %{
       P_d = normpdf(dst_ry,0,s_mj/3);
       
       dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
       
       idx = find(dst < 1);
       
       idx2 = find(dst_rypp > r/2 & dst < 1);
       
       P_bin12(i,idx) = single(P_d(idx));
       P_bin12(i,idx2) = single(normpdf(dst_rypp(idx2),r/2,s_mj/3));
       %}
       P_bin12(i,:) = single(P_d/sum(P_d));
   end
end

g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

% f = waitbar(0,'1','Name','Iterating through parameters for 8-16 Hz',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%{
gs_set = logspace(-7,1,100);
bs_set = logspace(-7,1,100);
err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,3);
gs_tmp = gs_set(idx);

gs_set = logspace(log10(min(gs_tmp)),log10(max(gs_tmp)),100);

[~,idx] = mink(err_bset,3);

bs_tmp = bs_set(idx);

bs_set = logspace(log10(min(bs_tmp)),log10(max(bs_tmp)),100);

err_gset = zeros(100,1);
err_bset = zeros(100,1);
for i=1:100
    g_tmp = sum((gs_set(i) * P_bin12),2);
    err_gset(i) =  (g_tmp-A{8})'*(g_tmp-A{8}); 
    b_tmp = sum((bs_set(i) * P_bin12),2);
    err_bset(i) =  (b_tmp-A{9})'*(b_tmp-A{9}); 
end

[~,idx] = mink(err_gset,1);

g_fi = single(gs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*gs_set(idx));

[~,idx] = mink(err_bset,1);

b_fi = single(bs_set(idx)*ones(n_grd,1) + rand(n_grd,1)*0.05*bs_set(idx));

% g_f = sum(bsxfun(@times,P_bin12,(A{8}.*A{6})))./sum(P_bin12);
% b_f = sum(bsxfun(@times,P_bin12,A{9}))./sum(P_bin12);

%%{

std_g = single(0.0025);
std_b = single(0.0025);

g_fmt = P_bin12 * g_fi;
b_fmt = P_bin12 * b_fi;

g_f = g_fi;
b_f = b_fi;

L_b = (sum(del2(b_f).^2)/n_grd)*n_a;
L_g = (sum(del2(g_f).^2)/n_grd)*n_a;

err_g = (g_fmt-A{8})'*(g_fmt-A{8}) + L_g;
err_b = (b_fmt-A{9})'*(b_fmt-A{9}) + L_b;

g_ff = single(0);
b_ff = single(0);

for iter=1:250000
   var = single(randi(n_grd));
   g_f2 = g_f;
   b_f2 = b_f;
   chus = randi(2);
   if chus == 1
       temp1 = g_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_g,(1-temp1)/std_g));
       g_f2(var) = temp1 + temp2 * std_g;
       L_g2 = (sum(del2(g_f2).^2)/n_grd)*n_a;
       g_fmt2 = P_bin12 * g_f2;
       err_g2 = (g_fmt2-A{8})'*(g_fmt2-A{8}) + L_g2;
       accp_prob = min(1,exp(-(err_g2-err_g)));
       if (accp_prob == 1 || accp_prob > rand(1))
           g_f = g_f2;
           err_g = err_g2;
       end
   else
       temp1 = b_f2(var);
       temp2 = single(trandn((0.001-temp1)/std_b,(1-temp1)/std_b));
       b_f2(var) = temp1 + temp2 * std_b;
       L_b2 = (sum(del2(b_f2).^2)/n_grd)*n_a;
       b_fmt2 = P_bin12 * b_f2;
       err_b2 = (b_fmt2-A{9})'*(b_fmt2-A{9}) + L_b2;
       accp_prob = min(1,exp(-(err_b2-err_b)));
       if (accp_prob == 1 || accp_prob > rand(1))
           b_f = b_f2;
           err_b = err_b2;
       end
   end
   if (iter > 200000)
      g_ff = (g_f + g_ff*(iter-200000-1))/(iter-200000);
      b_ff = (b_f + b_ff*(iter-200000-1))/(iter-200000);
      iter
   end
%    waitbar(iter/250000,f,sprintf('%f:%f:%f',iter,err_g,err_b));
%    if getappdata(f,'canceling')
%         break
%    end
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

g_f = g_ff;
b_f = b_ff;
%}

st1 = sprintf('%s/envnodes15g_16-32_%s.txt',wrk_path,sgmy);
st2 = sprintf('%s/envnodes15b_16-32_%s.txt',wrk_path,sgmy);
st3 = sprintf('%s/envnodes15a_16-32_%s.txt',wrk_path,sgmy);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');
fid5 = fopen(st3,'w');

Qsc_inv = g_f/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

albdo = Qsc_inv./(Qsc_inv + Qi_inv)*100;

un_arr = [bincord(:,1) bincord(:,2) log10(Qsc_inv')];
fprintf(fid3,'%f %f %f\n',un_arr');
fclose(fid3);
un_arr = [bincord(:,1) bincord(:,2) log10(Qi_inv')];
fprintf(fid4,'%f %f %f\n',un_arr');
fclose(fid4);
un_arr = [bincord(:,1) bincord(:,2) albdo'];
fprintf(fid5,'%f %f %f\n',un_arr');
fclose(fid5);