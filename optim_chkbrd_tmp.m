%% Assigning search directory for events and reading station coordinates and average station velocities 
clc;
close all;
clearvars;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

wkdir = '/home/user/env_EGEL/env_chktst_pkdbl_dfwln';

fldir = '/home/user/env_EGEL/env_gauss_7.2,3.6,1.8,.90_pkdbl';

mkdir(wkdir)

%% The grid dimensions for the mapping

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

fctr_mjr = 1/9; 
fctr_mnr = 1/9;

%% original inverted model

% abs_pdir = sprintf('%s/envnodes15g_1-2.txt',fldir);
% fid_rp = fopen(abs_pdir,'r');
% G_12 = textscan(fid_rp,'%f %f %f');
% fclose(fid_rp);
% 
% abs_pdir = sprintf('%s/envnodes15b_1-2.txt',fldir);
% fid_rp = fopen(abs_pdir,'r');
% B_12 = textscan(fid_rp,'%f %f %f');
% fclose(fid_rp);

%% For free anomaly test

%{
b_set = ones(n_grd,1)*0.01;
g_set = ones(n_grd,1)*0.0001;

xgv = [23.7891479, 23.7891479, 24.5891479, 24.5891479; ...% Milos
    25.1971436, 25.1971436, 25.9971436, 25.9971436; ...% Santorini
    26.9274902, 26.9274902, 27.7274902, 27.7274902; ...% Nisyros
    23.2429504, 23.2429504, 24.0429504, 24.0429504];% Methana

ygv = [37.2554499, 36.4554499, 36.4554499, 37.2554499; ...% Milos
    36.9261913, 36.1261913, 36.1261913, 36.9261913; ...% Santorini
    37.0915474, 36.2915474, 36.2915474, 37.0915474; ...% Nisyros
    38.1935076, 37.3935076, 37.3935076, 38.1935076];% Methana

for j = 1:length(ygv(:,1))
    polyg = [xgv(j,:)' ygv(j,:)'];
    [In1,On] = inpolygon(bincord(:,2),bincord(:,1),polyg(:,1),polyg(:,2));
    
    b_set(In1) = 0.1;
    g_set(In1) = 0.001;
end


chkm_dir = sprintf('%s/bin_chk_orig_bg.txt',wkdir);
fid = fopen(chkm_dir,'w');
fprintf(fid,'%f %f %f %f %f\n',[bincord log10(g_set) log10(b_set)]');
fclose(fid);
%}


%% for checkerboard test
%%{

Qi_invst = zeros(n_grd,1,'single');
Qs_invst = zeros(n_grd,1,'single');

xgv = 20.0:1.2:29.6;
ygv = 34.0:1.2:38.8;

for j = 1:length(ygv)-1
for i = 1:length(xgv)-1
polyg = [xgv(i),ygv(j); xgv(i+1),ygv(j); xgv(i+1),ygv(j+1); xgv(i),ygv(j+1)];
[In1,On] = inpolygon(bincord(:,2),bincord(:,1),polyg(:,1),polyg(:,2));
%In2 = inpolygon(bincord(pq2,2),bincord(pq2,1),polyg(:,1),polyg(:,2));
%In3 = inpolygon(bincord(pq3,2),bincord(pq3,1),polyg(:,1),polyg(:,2));
%In4 = inpolygon(bincord(pq4,2),bincord(pq4,1),polyg(:,1),polyg(:,2));
%In5 = inpolygon(bincord(pq5,2),bincord(pq5,1),polyg(:,1),polyg(:,2));
flag = (-1)^(i+j);
if flag > 0 
    Qi_invst(In1) = single(0.01);%10.^(B_12{3}(In1)) + 10.^(B_12{3}(In1))*0.05;
    Qs_invst(In1) = single(0.1);%10.^(G_12{3}(In1)) + 10.^(G_12{3}(In1))*0.05;
else
    Qi_invst(In1) = single(0.001);%10.^(B_12{3}(In1)) - 10.^(B_12{3}(In1))*0.05;
    Qs_invst(In1) = single(0.0001);%10.^(G_12{3}(In1)) - 10.^(G_12{3}(In1))*0.05;
end

end
end


chkm_dir = sprintf('%s/bin_chk_orig_bg.txt',wkdir);
fid = fopen(chkm_dir,'w');
%fprintf(fid,'%f %f %f %f %f\n',[bincord Qs_invst-10.^G_12{3} Qi_invst-10.^B_12{3}]');
fprintf(fid,'%f %f %f %f %f\n',[bincord Qs_invst Qi_invst]');
fclose(fid);

%}


%% considering 100 Hz sampling frequency

% Direct S-wave window length in samples
%d_S_len = single(1000);


%% Computing for 1-2 Hz band

abs_pdir = sprintf('%s/ray_prmtrs12.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
A = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/env_wn_1-2.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
B = textscan(fid_rp,'%f');
fclose(fid_rp);

wln_st = B{1};

% erth_rad = 6371;
% 
% arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);
% 
% r_n = sqrt(arclen.^2 + A{3}.^2);

b_orig = cell(length(A{1}),1);
g_orig = cell(length(A{1}),1);

b_inv = cell(length(A{1}),1);
g_st_inv = cell(length(A{1}),1);

P_bin12 = cell(length(A{1}),1);

freq = single(1.5);

g_tmp = ((Qs_invst))*(2*pi*freq);
b_tmp = ((Qi_invst))*(2*pi*freq);

f = waitbar(0,'1','Name','Iterating through files for 1-2 Hz',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

g_rng = single(logspace(-7,1,100));

tic

rng('shuffle');

for i=1:length(A{1})
    ev_lat = single(A{1}(i));
    ev_lon = single(A{2}(i));
    ev_dep = single(A{3}(i));
    
    v_avg = single(A{6}(i));
    t_end = single(A{7}(i));
    sz = round(t_end * 100);
    
    st_lat = single(A{4}(i));
    st_lon = single(A{5}(i));
    
    [ev_x,ev_y,ev_z] = geodetic2enu(ev_lat,ev_lon,0,34,18,0,wgs84);
    
    W = single(rand*exp(17));
    
    R = single(rand*10);
    
    [st_x,st_y,st_z] = geodetic2enu(st_lat,st_lon,0,34,18,0,wgs84);
    r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
    if (1)%(r < 60)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s_mj = (v_avg * t_end + r)/2;
        s_mn = sqrt(4*s_mj^2 - r^2)/2;
        Cx = (st_x + ev_x)/2;
        Cy = (st_y + ev_y)/2;
        slp_ry = (st_y-ev_y)/(st_x-ev_x);
        dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
        slp_rypp = -1/slp_ry;
        dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
        
        P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
        P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
        
%         dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
%         
%         idx = dst > 1;
%         
%         P_d(idx) = 0;
%         idx = find(P_d~=0);
        
        tempg = (g_tmp .* P_d)./sum(P_d);
        tempb = (b_tmp .* P_d)./sum(P_d);
        
        g_st = single(sum(tempg)/v_avg);
        b_vl = single(sum(tempb));
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        temp2 = g_func_env_tst(g_st,v_avg,r,sz,...
            b_vl);
        
        temp2 = temp2.*10^(-9);
        
        %%{
        for k = 1:length(temp2)
            temp2(k) = temp2(k) + ((-1).^(round(rand*100)))* temp2(k) *0.10 * rand;
        end
        %}
        %figure(1); plot(temp2(2:end));
%         idx_tmp = single(find(temp2(1:end-800)==max(temp2(1:end-800))));
%         S_idx = idx_tmp(1);
%         
%         w_len = S_idx*2;
        %S_idx= 500;
        
        w_len = wln_st(i);
        
        S_idx = round(w_len/2);
        
        dcrs_flg = 1;
        
%         if (S_idx > 10)
% 
%             w_len = 500;
%             S_idx = 250;
%             dcrs_flg = 0;
%         end
        
        E = zeros(sz-w_len,1,'single');
        G_f = zeros(sz-w_len,1,'single');
    
    
        E_f = temp2*R*W;
        
        
        tm = [single((S_idx)/100) ...
            single(linspace(((w_len)/100+0.01),t_end,sz-w_len-1))];
        
        E(1) = mean(E_f(1:w_len+1));
        E(2:(sz-w_len)) = E_f((w_len+2):end);
        
        E = log(E);
        
        %%{
        err_sum = g_func_wrpr_mtrx(g_rng,w_len,E,...
            tm,v_avg,r,sz);
        
        [~,g_idx] = mink(err_sum,1);
        
        if (g_idx < 100) && (g_idx > 1)
            
            g_rng2 = logspace(log10(g_rng(g_idx-1)),...
                log10(g_rng(g_idx+1)),100);
            
            err_sum = g_func_wrpr_mtrx(g_rng2,w_len,E,...
                tm,v_avg,r,sz);
            
            [~,g_idx2] = mink(err_sum,1);
            
            g_st_w = g_rng2(g_idx2);
            
        else
            
            g_st_w = g_rng(g_idx);
        end
    
        %}
        
        temp3 = g_func_env_tst(g_st_w,v_avg,r,sz,...
            0);
        
        temp1 = temp3.*10^(-9);
        
        G_f(1) = mean(temp1(1:w_len+1));
        G_f(2:(sz-w_len)) = temp1((w_len+2):end);
        
        G_f = log(G_f);
        
        X = [ones((sz-w_len),1,'single') -tm'];
        
        b = X\(E-G_f);
        
        E_calc = G_f - b(2).*tm' + b(1);
        
        sum_tmp = sum(abs(exp(E)-exp(E_calc)))/length(E) * 100;
        err_tmp = abs(exp(E(1)) - exp(E_calc(1)));
        error_prom = abs(err_tmp/exp(E(1)))*100;
        
        if(b(2) > 0) && (dcrs_flg) && (error_prom < 1)
            
            b_inv{i} = single(b(2)/(2*pi*freq));
            g_st_inv{i} = single(g_st_w*v_avg/(2*pi*freq));
            
            P_bin12{i} = single(P_d'/sum(P_d));
            g_orig{i} = single(g_st);
            b_orig{i} = single(b_vl);
            
            %{
            i
            g_st
            g_st_w
            b_vl
            b(2)
            log(W*R)
            b(1)
            error_prom
            log10(sum_tmp)
        %}
        else
            b_inv{i} = single(b_inv{i});
            g_st_inv{i} = single(g_st_inv{i});
            
            P_bin12{i} = single(P_bin12{i});
            g_orig{i} = single(g_orig{i});
            b_orig{i} = single(b_orig{i});
        end
%         plot(E_calc); hold on; plot(E); hold off;
        
        if (rem(i,100)==0)
        waitbar(i/length(A{1}),f,sprintf('%12.0f',i))
        end
        %{
        i
        w_len
        g_st
        g_st_w
        b_vl
        b(2)
        log(W*R)
        b(1)
        error_prom
        sum_tmp
        %}
        
        %{
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                t_fnl = linspace(0,tm(end),length(E_f));
                
                h1 = figure(i*15); set(gca,'FontName','Helvetica',...
                    'Fontsize',14,'FontWeight','bold'); hold on;
                scatter(tm(1),E(1),'r','LineWidth',2); ...
                    plot(tm,E,'r','LineWidth',2); %legend(string(r));
                hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2);...
                    plot(tm,E_calc,'b','LineWidth',2); xlabel('Time (in s)');...
                    ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}',...
                    'S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                
                figure(i*150); set(gca,'FontName','Helvetica','Fontsize',...
                    14,'FontWeight','bold'); hold on;...
                    xlabel('Time (in s)');...
                    ylabel('Amplitude');
                plot(t_fnl,E_f,'LineWidth',2); legend('Smoothed S-wave env.');...
                    hold off;
        %}
    else
        b_inv{i} = single(b_inv{i});
        g_st_inv{i} = single(g_st_inv{i});
        
        P_bin12{i} = single(P_bin12{i});
        g_orig{i} = single(g_orig{i});
        b_orig{i} = single(b_orig{i});
    end
    
end

P_bin12 = cell2mat(P_bin12);
g_st_inv = cell2mat(g_st_inv);
b_inv = cell2mat(b_inv);
g_orig = cell2mat(g_orig);
b_orig = cell2mat(b_orig);

toc

% delete(f);

%{
g_f = lsqlin(double(P_bin12),double((g_st_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);
b_f = lsqlin(double(P_bin12),double((b_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);

 %%%%%%%%%%%%%%%%%%%%%% Or %%%%%%%%%%%%%%%%%%%
%}

%{

g_f = P_bin12\g_st_inv;

b_f = P_bin12\b_inv;

idx = find(g_f > 0 & b_f > 0);

arset = [bincord(idx,:) log10(g_f(idx)) log10(b_f(idx))];

%}

%%{
g_f = sum(bsxfun(@times,g_st_inv,P_bin12))...
    ./sum(P_bin12);
b_f = sum(bsxfun(@times,b_inv,P_bin12))...
    ./sum(P_bin12);

%arset = [bincord g_f'-10.^G_12{3} b_f'-10.^B_12{3}];
arset = [bincord g_f' b_f'];
%}

save('vars_orig_chkbrd12.mat','b_orig','g_orig');
save('vars_inv_chkbrd12.mat','b_inv','g_st_inv');
save('prob_chkbrd12.mat','P_bin12');

sp = sprintf('%s/bin_chk_31_46_1_1-2Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f\n',arset');
fclose(fid);
%}
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

diffb = abs(b_orig - b_inv)./b_orig * 100;
diffg = abs(g_orig - g_st_inv)./g_orig * 100;
figure(1); histogram(diffg);
figure(2); histogram(diffb);

%%{

% clear b_orig g_orig lgW_orig R_orig ...
%   b_inv g_st_inv lgW_inv R_inv P_bin12

%% Computing for 2-4 Hz band

abs_pdir = sprintf('%s/ray_prmtrs24.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
A = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/env_wn_2-4.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
B = textscan(fid_rp,'%f');
fclose(fid_rp);

wln_st = B{1};


% erth_rad = 6371;
% 
% arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);
% 
% r_n = sqrt(arclen.^2 + A{3}.^2);

b_orig = cell(length(A{1}),1);
g_orig = cell(length(A{1}),1);

b_inv = cell(length(A{1}),1);
g_st_inv = cell(length(A{1}),1);

P_bin12 = cell(length(A{1}),1);

freq = single(3);

g_tmp = ((Qs_invst));%*(2*pi*freq);
b_tmp = ((Qi_invst));%*(2*pi*freq);

f = waitbar(0,'1','Name','Iterating through files for 2-4 Hz',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

g_rng = single(logspace(-7,1,100));

tic

rng('shuffle');

for i=1:length(A{1})
    ev_lat = single(A{1}(i));
    ev_lon = single(A{2}(i));
    ev_dep = single(A{3}(i));
    
    v_avg = single(A{6}(i));
    t_end = single(A{7}(i));
    sz = round(t_end * 100);
    
    st_lat = single(A{4}(i));
    st_lon = single(A{5}(i));
    
    [ev_x,ev_y,ev_z] = geodetic2enu(ev_lat,ev_lon,0,34,18,0,wgs84);
    
    W = single(rand*exp(17));
    
    R = single(rand*10);
    
    [st_x,st_y,st_z] = geodetic2enu(st_lat,st_lon,0,34,18,0,wgs84);
    r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
    if (1)%(r < 60)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s_mj = (v_avg * t_end + r)/2;
        s_mn = sqrt(4*s_mj^2 - r^2)/2;
        Cx = (st_x + ev_x)/2;
        Cy = (st_y + ev_y)/2;
        slp_ry = (st_y-ev_y)/(st_x-ev_x);
        dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
        slp_rypp = -1/slp_ry;
        dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
        
        P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
        P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
        
%         dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
%         
%         idx = dst > 1;
%         
%         P_d(idx) = 0;
%         idx = find(P_d~=0);
        
        tempg = (g_tmp .* P_d)./sum(P_d);
        tempb = (b_tmp .* P_d)./sum(P_d);
        
        g_st = single(sum(tempg)/v_avg);
        b_vl = single(sum(tempb));
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        temp2 = g_func_env_tst(g_st,v_avg,r,sz,...
            b_vl);
        
        temp2 = temp2.*10^(-9);
        
        %%{
        for k = 1:length(temp2)
            temp2(k) = temp2(k) + ((-1).^(round(rand*100)))* temp2(k) *0.10 * rand;
        end
        %}
        %figure(1); plot(temp2(2:end));
        
%         idx_tmp = single(find(temp2(1:end-800)==max(temp2(1:end-800))));
%         S_idx = idx_tmp(1);
%         
%         w_len = S_idx*2;

        w_len = wln_st(i);
        S_idx = round(w_len/2);
        
        dcrs_flg = 1;
        
%         if (S_idx > 10)
% 
%             w_len = 500;
%             S_idx = 250;
%             dcrs_flg = 0;
%         end
        
        E = zeros(sz-w_len,1,'single');
        G_f = zeros(sz-w_len,1,'single');
    
    
        E_f = temp2*R*W;
        
        tm = [single((S_idx)/100) ...
            single(linspace(((w_len)/100+0.01),t_end,sz-w_len-1))];
        
        E(1) = mean(E_f(1:w_len+1));
        E(2:(sz-w_len)) = E_f((w_len+2):end);
        
        E = log(E);
        
        %%{
        err_sum = g_func_wrpr_mtrx(g_rng,w_len,E,...
            tm,v_avg,r,sz);
        
        [~,g_idx] = mink(err_sum,1);
        
        if (g_idx < 100) && (g_idx > 1)
            
            g_rng2 = logspace(log10(g_rng(g_idx-1)),...
                log10(g_rng(g_idx+1)),100);
            
            err_sum = g_func_wrpr_mtrx(g_rng2,w_len,E,...
                tm,v_avg,r,sz);
            
            [~,g_idx2] = mink(err_sum,1);
            
            g_st_w = g_rng2(g_idx2);
            
        else
            
            g_st_w = g_rng(g_idx);
        end
    
        %}
        
        temp3 = g_func_env_tst(g_st_w,v_avg,r,sz,...
            0);
        
        temp1 = temp3.*10^(-9);
        
        G_f(1) = mean(temp1(1:w_len+1));
        G_f(2:(sz-w_len)) = temp1((w_len+2):end);
        
        G_f = log(G_f);
        
        X = [ones((sz-w_len),1,'single') -tm'];
        
        b = X\(E-G_f);
        
        E_calc = G_f - b(2).*tm' + b(1);
        
        sum_tmp = sum(abs((E(2:end)-E_calc(2:end))))/length(E) * 100;
        err_tmp = abs(exp(E(1)) - exp(E_calc(1)));
        error_prom = abs(err_tmp/exp(E(1)))*100;
        
        if(b(2) > 0) && (dcrs_flg) && (error_prom < 1)
            
            b_inv{i} = single(b(2));%/(2*pi*freq));
            g_st_inv{i} = single(g_st_w*v_avg);%/(2*pi*freq));
            
            P_bin12{i} = single(P_d'/sum(P_d));
            g_orig{i} = single(g_st);
            b_orig{i} = single(b_vl);
        else
            b_inv{i} = single(b_inv{i});
            g_st_inv{i} = single(g_st_inv{i});
            
            P_bin12{i} = single(P_bin12{i});
            g_orig{i} = single(g_orig{i});
            b_orig{i} = single(b_orig{i});
        end
        %plot(E_calc); hold on; plot(E); hold off;
        
        if (rem(i,100)==0)
        waitbar(i/length(A{1}),f,sprintf('%12.0f',i))
        end
        %{
        i
        w_len
        g_st
        g_st_w
        b_vl
        b(2)
        log(W*R)
        b(1)
        error_prom
        %}
        
        %{
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                t_fnl = linspace(-1.0,tm(end),sz-1);
                
                h1 = figure(j*15); set(gca,'FontName','Helvetica',...
                    'Fontsize',14,'FontWeight','bold'); hold on;
                scatter(tm(1),E(1),'r','LineWidth',2); ...
                    plot(tm,E,'r','LineWidth',2); %legend(string(r));
                hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2);...
                    plot(tm,E_calc,'b','LineWidth',2); ...
                    plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                    ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}',...
                    'S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                
                figure(j*150); set(gca,'FontName','Helvetica','Fontsize',...
                    14,'FontWeight','bold'); hold on;...
                    xlabel('Time (in s)');...
                    ylabel('Amplitude');
                plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.');...
                    hold off;
        %}
    else
        b_inv{i} = single(b_inv{i});
        g_st_inv{i} = single(g_st_inv{i});
        
        P_bin12{i} = single(P_bin12{i});
        g_orig{i} = single(g_orig{i});
        b_orig{i} = single(b_orig{i});
    end
    
end

P_bin12 = cell2mat(P_bin12);
g_st_inv = cell2mat(g_st_inv);
b_inv = cell2mat(b_inv);
g_orig = cell2mat(g_orig);
b_orig = cell2mat(b_orig);

toc

% delete(f);

%{
g_f = lsqlin(double(P_bin12),double((g_st_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);
b_f = lsqlin(double(P_bin12),double((b_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);

 %%%%%%%%%%%%%%%%%%%%%% Or %%%%%%%%%%%%%%%%%%%
%}

%{

g_f = P_bin12\g_st_inv;

b_f = P_bin12\b_inv;

idx = find(g_f > 0 & b_f > 0);

arset = [bincord(idx,:) log10(g_f(idx)) log10(b_f(idx))];

%}

%%{
g_f = sum(bsxfun(@times,g_st_inv,P_bin12))...
    ./sum(P_bin12);
b_f = sum(bsxfun(@times,b_inv,P_bin12))...
    ./sum(P_bin12);

%arset = [bincord g_f'-10.^G_12{3} b_f'-10.^B_12{3}];
arset = [bincord g_f' b_f'];
%}

save('vars_orig_chkbrd24.mat','b_orig','g_orig');
save('vars_inv_chkbrd24.mat','b_inv','g_st_inv');
save('prob_chkbrd24.mat','P_bin12');

sp = sprintf('%s/bin_chk_31_46_1_2-4Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f\n',arset');
fclose(fid);
%}
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

diffb = abs(b_orig - b_inv)./b_orig * 100;
diffg = abs(g_orig - g_st_inv)./g_orig * 100;
figure(1); histogram(diffg);
figure(2); histogram(diffb);

%% Computing for 4-8 Hz band

abs_pdir = sprintf('%s/ray_prmtrs48.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
A = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/env_wn_4-8.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
B = textscan(fid_rp,'%f');
fclose(fid_rp);

wln_st = B{1};


% erth_rad = 6371;
% 
% arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);
% 
% r_n = sqrt(arclen.^2 + A{3}.^2);

b_orig = cell(length(A{1}),1);
g_orig = cell(length(A{1}),1);

b_inv = cell(length(A{1}),1);
g_st_inv = cell(length(A{1}),1);

P_bin12 = cell(length(A{1}),1);

freq = single(6);

g_tmp = ((Qs_invst));%*(2*pi*freq);
b_tmp = ((Qi_invst));%*(2*pi*freq);

f = waitbar(0,'1','Name','Iterating through files for 4-8 Hz',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

g_rng = single(logspace(-7,1,100));

tic

rng('shuffle');

for i=1:length(A{1})
    ev_lat = single(A{1}(i));
    ev_lon = single(A{2}(i));
    ev_dep = single(A{3}(i));
    
    v_avg = single(A{6}(i));
    t_end = single(A{7}(i));
    sz = round(t_end * 100);
    
    st_lat = single(A{4}(i));
    st_lon = single(A{5}(i));
    
    [ev_x,ev_y,ev_z] = geodetic2enu(ev_lat,ev_lon,0,34,18,0,wgs84);
    
    W = single(rand*exp(17));
    
    R = single(rand*10);
    
    [st_x,st_y,st_z] = geodetic2enu(st_lat,st_lon,0,34,18,0,wgs84);
    r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
    if (1)%(r < 60)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s_mj = (v_avg * t_end + r)/2;
        s_mn = sqrt(4*s_mj^2 - r^2)/2;
        Cx = (st_x + ev_x)/2;
        Cy = (st_y + ev_y)/2;
        slp_ry = (st_y-ev_y)/(st_x-ev_x);
        dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
        slp_rypp = -1/slp_ry;
        dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
        
        P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
        P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
        
%         dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
%         
%         idx = dst > 1;
%         
%         P_d(idx) = 0;
%         idx = find(P_d~=0);
        
        tempg = (g_tmp .* P_d)./sum(P_d);
        tempb = (b_tmp .* P_d)./sum(P_d);
        
        g_st = single(sum(tempg)/v_avg);
        b_vl = single(sum(tempb));
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        temp2 = g_func_env_tst(g_st,v_avg,r,sz,...
            b_vl);
        
        temp2 = temp2.*10^(-9);
        
        %%{
        for k = 1:length(temp2)
            temp2(k) = temp2(k) + ((-1).^(round(rand*100)))* temp2(k) *0.10 * rand;
        end
        %}
        %figure(1); plot(temp2(2:end));
        
        %         idx_tmp = single(find(temp2(1:end-800)==max(temp2(1:end-800))));
%         S_idx = idx_tmp(1);
%         
%         w_len = S_idx*2;

        w_len = wln_st(i);
        S_idx = round(w_len/2);
        
        dcrs_flg = 1;
        
%         if (S_idx > 10)
% 
%             w_len = 500;
%             S_idx = 250;
%             dcrs_flg = 0;
%         end
        
        E = zeros(sz-w_len,1,'single');
        G_f = zeros(sz-w_len,1,'single');
    
    
        E_f = temp2*R*W;
        
        tm = [single((S_idx)/100) ...
            single(linspace(((w_len)/100+0.01),t_end,sz-w_len-1))];
        
        E(1) = mean(E_f(1:w_len+1));
        E(2:(sz-w_len)) = E_f((w_len+2):end);
        
        E = log(E);
        
        %%{
        err_sum = g_func_wrpr_mtrx(g_rng,w_len,E,...
            tm,v_avg,r,sz);
        
        [~,g_idx] = mink(err_sum,1);
        
        if (g_idx < 100) && (g_idx > 1)
            
            g_rng2 = logspace(log10(g_rng(g_idx-1)),...
                log10(g_rng(g_idx+1)),100);
            
            err_sum = g_func_wrpr_mtrx(g_rng2,w_len,E,...
                tm,v_avg,r,sz);
            
            [~,g_idx2] = mink(err_sum,1);
            
            g_st_w = g_rng2(g_idx2);
            
        else
            
            g_st_w = g_rng(g_idx);
        end
    
        %}
        
        temp3 = g_func_env_tst(g_st_w,v_avg,r,sz,...
            0);
        
        temp1 = temp3.*10^(-9);
        
        G_f(1) = mean(temp1(1:w_len+1));
        G_f(2:(sz-w_len)) = temp1((w_len+2):end);
        
        G_f = log(G_f);
        
        X = [ones((sz-w_len),1,'single') -tm'];
        
        b = X\(E-G_f);
        
        E_calc = G_f - b(2).*tm' + b(1);
        
        sum_tmp = sum(abs((E(2:end)-E_calc(2:end))))/length(E) * 100;
        err_tmp = abs(exp(E(1)) - exp(E_calc(1)));
        error_prom = abs(err_tmp/exp(E(1)))*100;
        
        if(b(2) > 0) && (dcrs_flg) && (error_prom < 1)
            
            b_inv{i} = single(b(2));%/(2*pi*freq));
            g_st_inv{i} = single(g_st_w*v_avg);%/(2*pi*freq));
            
            P_bin12{i} = single(P_d'/sum(P_d));
            g_orig{i} = single(g_st);
            b_orig{i} = single(b_vl);
        else
            b_inv{i} = single(b_inv{i});
            g_st_inv{i} = single(g_st_inv{i});
            
            P_bin12{i} = single(P_bin12{i});
            g_orig{i} = single(g_orig{i});
            b_orig{i} = single(b_orig{i});
        end
        %plot(E_calc); hold on; plot(E); hold off;
        
        if (rem(i,100)==0)
        waitbar(i/length(A{1}),f,sprintf('%12.0f',i))
        end
        %{
        i
        w_len
        g_st
        g_st_w
        b_vl
        b(2)
        log(W*R)
        b(1)
        error_prom
        %}
        
        %{
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                t_fnl = linspace(-1.0,tm(end),sz-1);
                
                h1 = figure(j*15); set(gca,'FontName','Helvetica',...
                    'Fontsize',14,'FontWeight','bold'); hold on;
                scatter(tm(1),E(1),'r','LineWidth',2); ...
                    plot(tm,E,'r','LineWidth',2); %legend(string(r));
                hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2);...
                    plot(tm,E_calc,'b','LineWidth',2); ...
                    %plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                    ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}',...
                    'S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                
                figure(j*150); set(gca,'FontName','Helvetica','Fontsize',...
                    14,'FontWeight','bold'); hold on;...
                    xlabel('Time (in s)');...
                    ylabel('Amplitude');
                plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.');...
                    hold off;
        %}
    else
        b_inv{i} = single(b_inv{i});
        g_st_inv{i} = single(g_st_inv{i});
        
        P_bin12{i} = single(P_bin12{i});
        g_orig{i} = single(g_orig{i});
        b_orig{i} = single(b_orig{i});
    end
    
end

P_bin12 = cell2mat(P_bin12);
g_st_inv = cell2mat(g_st_inv);
b_inv = cell2mat(b_inv);
g_orig = cell2mat(g_orig);
b_orig = cell2mat(b_orig);

toc

% delete(f);

%{
g_f = lsqlin(double(P_bin12),double((g_st_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);
b_f = lsqlin(double(P_bin12),double((b_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);

 %%%%%%%%%%%%%%%%%%%%%% Or %%%%%%%%%%%%%%%%%%%
%}

%{

g_f = P_bin12\g_st_inv;

b_f = P_bin12\b_inv;

idx = find(g_f > 0 & b_f > 0);

arset = [bincord(idx,:) log10(g_f(idx)) log10(b_f(idx))];

%}

%%{
g_f = sum(bsxfun(@times,g_st_inv,P_bin12))...
    ./sum(P_bin12);
b_f = sum(bsxfun(@times,b_inv,P_bin12))...
    ./sum(P_bin12);

%arset = [bincord g_f'-10.^G_12{3} b_f'-10.^B_12{3}];
arset = [bincord g_f' b_f'];
%}

save('vars_orig_chkbrd48.mat','b_orig','g_orig');
save('vars_inv_chkbrd48.mat','b_inv','g_st_inv');
save('prob_chkbrd48.mat','P_bin12');

sp = sprintf('%s/bin_chk_31_46_1_4-8Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f\n',arset');
fclose(fid);
%}
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

diffb = abs(b_orig - b_inv)./b_orig * 100;
diffg = abs(g_orig - g_st_inv)./g_orig * 100;
figure(1); histogram(diffg);
figure(2); histogram(diffb);

%% Computing for 8-16 Hz band

abs_pdir = sprintf('%s/ray_prmtrs816.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
A = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
fclose(fid_rp);

abs_pdir = sprintf('%s/env_wn_8-16.txt',fldir);
fid_rp = fopen(abs_pdir,'r');
B = textscan(fid_rp,'%f');
fclose(fid_rp);

wln_st = B{1};


% erth_rad = 6371;
% 
% arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);
% 
% r_n = sqrt(arclen.^2 + A{3}.^2);

b_orig = cell(length(A{1}),1);
g_orig = cell(length(A{1}),1);

b_inv = cell(length(A{1}),1);
g_st_inv = cell(length(A{1}),1);

P_bin12 = cell(length(A{1}),1);

freq = single(12);

g_tmp = ((Qs_invst));%*(2*pi*freq);
b_tmp = ((Qi_invst));%*(2*pi*freq);

f = waitbar(0,'1','Name','Iterating through files for 8-16 Hz',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

g_rng = single(logspace(-7,1,100));

tic

rng('shuffle');

for i=1:length(A{1})
    ev_lat = single(A{1}(i));
    ev_lon = single(A{2}(i));
    ev_dep = single(A{3}(i));
    
    v_avg = single(A{6}(i));
    t_end = single(A{7}(i));
    sz = round(t_end * 100);
    
    st_lat = single(A{4}(i));
    st_lon = single(A{5}(i));
    
    [ev_x,ev_y,ev_z] = geodetic2enu(ev_lat,ev_lon,0,34,18,0,wgs84);
    
    W = single(rand*exp(17));
    
    R = single(rand*10);
    
    [st_x,st_y,st_z] = geodetic2enu(st_lat,st_lon,0,34,18,0,wgs84);
    r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
    if (1)%(r < 60)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s_mj = (v_avg * t_end + r)/2;
        s_mn = sqrt(4*s_mj^2 - r^2)/2;
        Cx = (st_x + ev_x)/2;
        Cy = (st_y + ev_y)/2;
        slp_ry = (st_y-ev_y)/(st_x-ev_x);
        dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
        slp_rypp = -1/slp_ry;
        dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
        
        P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
        P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
        
%         dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
%         
%         idx = dst > 1;
%         
%         P_d(idx) = 0;
%         idx = find(P_d~=0);
        
        tempg = (g_tmp .* P_d)./sum(P_d);
        tempb = (b_tmp .* P_d)./sum(P_d);
        
        g_st = single(sum(tempg)/v_avg);
        b_vl = single(sum(tempb));
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        temp2 = g_func_env_tst(g_st,v_avg,r,sz,...
            b_vl);
        
        temp2 = temp2.*10^(-9);
        
        %%{
        for k = 1:length(temp2)
            temp2(k) = temp2(k) + ((-1).^(round(rand*100)))* temp2(k) *0.10 * rand;
        end
        %}
        %figure(1); plot(temp2(2:end));
        
        %         idx_tmp = single(find(temp2(1:end-800)==max(temp2(1:end-800))));
%         S_idx = idx_tmp(1);
%         
%         w_len = S_idx*2;

        w_len = wln_st(i);
        S_idx = round(w_len/2);
        
        dcrs_flg = 1;
        
%         if (S_idx > 10)
% 
%             w_len = 500;
%             S_idx = 250;
%             dcrs_flg = 0;
%         end
        
        E = zeros(sz-w_len,1,'single');
        G_f = zeros(sz-w_len,1,'single');
    
    
        E_f = temp2*R*W;
        
        tm = [single((S_idx)/100) ...
            single(linspace(((w_len)/100+0.01),t_end,sz-w_len-1))];
        
        E(1) = mean(E_f(1:w_len+1));
        E(2:(sz-w_len)) = E_f((w_len+2):end);
        
        E = log(E);
        
        %%{
        err_sum = g_func_wrpr_mtrx(g_rng,w_len,E,...
            tm,v_avg,r,sz);
        
        [~,g_idx] = mink(err_sum,1);
        
        if (g_idx < 100) && (g_idx > 1)
            
            g_rng2 = logspace(log10(g_rng(g_idx-1)),...
                log10(g_rng(g_idx+1)),100);
            
            err_sum = g_func_wrpr_mtrx(g_rng2,w_len,E,...
                tm,v_avg,r,sz);
            
            [~,g_idx2] = mink(err_sum,1);
            
            g_st_w = g_rng2(g_idx2);
            
        else
            
            g_st_w = g_rng(g_idx);
        end
    
        %}
        
        temp3 = g_func_env_tst(g_st_w,v_avg,r,sz,...
            0);
        
        temp1 = temp3.*10^(-9);
        
        G_f(1) = mean(temp1(1:w_len+1));
        G_f(2:(sz-w_len)) = temp1((w_len+2):end);
        
        G_f = log(G_f);
        
        X = [ones((sz-w_len),1,'single') -tm'];
        
        b = X\(E-G_f);
        
        E_calc = G_f - b(2).*tm' + b(1);
        
        sum_tmp = sum(abs((E(2:end)-E_calc(2:end))))/length(E) * 100;
        err_tmp = abs(exp(E(1)) - exp(E_calc(1)));
        error_prom = abs(err_tmp/exp(E(1)))*100;
        
        if(b(2) > 0) && (dcrs_flg) && (error_prom < 1)
            
            b_inv{i} = single(b(2));%/(2*pi*freq));
            g_st_inv{i} = single(g_st_w*v_avg);%/(2*pi*freq));
            
            P_bin12{i} = single(P_d'/sum(P_d));
            g_orig{i} = single(g_st);
            b_orig{i} = single(b_vl);
            
            %{
        g_st
        g_st_w
        b_vl
        b(2)
        log(W*R)
        b(1)
        error_prom
        %}
        
        else
            b_inv{i} = single(b_inv{i});
            g_st_inv{i} = single(g_st_inv{i});
            
            P_bin12{i} = single(P_bin12{i});
            g_orig{i} = single(g_orig{i});
            b_orig{i} = single(b_orig{i});
        end
        %plot(E_calc); hold on; plot(E); hold off;
        
        if (rem(i,100)==0)
        waitbar(i/length(A{1}),f,sprintf('%12.0f',i))
        end
        %{
        i
        w_len
        g_st
        g_st_w
        b_vl
        b(2)
        log(W*R)
        b(1)
        error_prom
        %}
        
        %{
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                t_fnl = linspace(-1.0,tm(end),sz-1);
                
                h1 = figure(j*15); set(gca,'FontName','Helvetica',...
                    'Fontsize',14,'FontWeight','bold'); hold on;
                scatter(tm(1),E(1),'r','LineWidth',2); ...
                    plot(tm,E,'r','LineWidth',2); %legend(string(r));
                hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2);...
                    plot(tm,E_calc,'b','LineWidth',2); ...
                    plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                    ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}',...
                    'S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                
                figure(j*150); set(gca,'FontName','Helvetica','Fontsize',...
                    14,'FontWeight','bold'); hold on;...
                    xlabel('Time (in s)');...
                    ylabel('Amplitude');
                plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.');...
                    hold off;
        %}
    else
        b_inv{i} = single(b_inv{i});
        g_st_inv{i} = single(g_st_inv{i});
        
        P_bin12{i} = single(P_bin12{i});
        g_orig{i} = single(g_orig{i});
        b_orig{i} = single(b_orig{i});
    end
    
end

P_bin12 = cell2mat(P_bin12);
g_st_inv = cell2mat(g_st_inv);
b_inv = cell2mat(b_inv);
g_orig = cell2mat(g_orig);
b_orig = cell2mat(b_orig);

toc

% delete(f);

%{
g_f = lsqlin(double(P_bin12),double((g_st_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);
b_f = lsqlin(double(P_bin12),double((b_inv)),[],[],[],[],...
    ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);

 %%%%%%%%%%%%%%%%%%%%%% Or %%%%%%%%%%%%%%%%%%%
%}

%{

g_f = P_bin12\g_st_inv;

b_f = P_bin12\b_inv;

idx = find(g_f > 0 & b_f > 0);

arset = [bincord(idx,:) log10(g_f(idx)) log10(b_f(idx))];

%}

%%{
g_f = sum(bsxfun(@times,g_st_inv,P_bin12))...
    ./sum(P_bin12);
b_f = sum(bsxfun(@times,b_inv,P_bin12))...
    ./sum(P_bin12);

%arset = [bincord g_f'-10.^G_12{3} b_f'-10.^B_12{3}];
arset = [bincord g_f' b_f'];
%}

save('vars_orig_chkbrd816.mat','b_orig','g_orig');
save('vars_inv_chkbrd816.mat','b_inv','g_st_inv');
save('prob_chkbrd816.mat','P_bin12');

sp = sprintf('%s/bin_chk_31_46_1_8-16Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f\n',arset');
fclose(fid);
%}
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

diffb = abs(b_orig - b_inv)./b_orig * 100;
diffg = abs(g_orig - g_st_inv)./g_orig * 100;
figure(1); histogram(diffg);
figure(2); histogram(diffb);

%% Computing for 16-32 Hz band

% abs_pdir = sprintf('%s/ray_prmtrs1632.txt',fldir);
% fid_rp = fopen(abs_pdir,'r');
% A = textscan(fid_rp,'%f %f %f %f %f %f %f %f %f');
% fclose(fid_rp);
% 
% % erth_rad = 6371;
% % 
% % arclen = distance([A{1} A{2}],[A{4} A{5}],erth_rad);
% % 
% % r_n = sqrt(arclen.^2 + A{3}.^2);
% 
% b_orig = cell(length(A{1}),1);
% g_orig = cell(length(A{1}),1);
% 
% b_inv = cell(length(A{1}),1);
% g_st_inv = cell(length(A{1}),1);
% 
% P_bin12 = cell(length(A{1}),1);
% 
% freq = single(24);
% 
% g_tmp = ((Qs_invst))*(2*pi*freq);
% b_tmp = ((Qi_invst))*(2*pi*freq);
% 
% f = waitbar(0,'1','Name','Iterating through files for 16-32 Hz',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% 
% g_rng = single(logspace(-7,1,100));
% 
% tic
% 
% rng('shuffle');
% 
% for i=1:length(A{1})
%     ev_lat = A{1}(i);
%     ev_lon = A{2}(i);
%     ev_dep = A{3}(i);
%     
%     v_avg = A{6}(i);
%     t_end = A{7}(i);
%     sz = round(t_end * 100);
%     
%     st_lat = A{4}(i);
%     st_lon = A{5}(i);
%     
%     [ev_x,ev_y,ev_z] = geodetic2enu(ev_lat,ev_lon,0,34,18,0,wgs84);
%     
%     W = single(rand*exp(17));
%     
%     R = single(rand*10);
%     
%     [st_x,st_y,st_z] = geodetic2enu(st_lat,st_lon,0,34,18,0,wgs84);
%     r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
%     if (1)%(r < 60)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         s_mj = (v_avg * t_end + r)/2;
%         s_mn = sqrt(4*s_mj^2 - r^2)/2;
%         Cx = (st_x + ev_x)/2;
%         Cy = (st_y + ev_y)/2;
%         slp_ry = (st_y-ev_y)/(st_x-ev_x);
%         dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
%         slp_rypp = -1/slp_ry;
%         dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
%         
%         P_sgma = [s_mj^2*fctr_mjr s_mn^2*fctr_mnr];
%         P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
%         
% %         dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
% %         
% %         idx = dst > 1;
% %         
% %         P_d(idx) = 0;
% %         idx = find(P_d~=0);
%         
%         tempg = (g_tmp .* P_d)./sum(P_d);
%         tempb = (b_tmp .* P_d)./sum(P_d);
%         
%         g_st = sum(tempg)/v_avg;
%         b_vl = sum(tempb);
%     
%     
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         w_len = single(d_S_len);
%         
%         E = zeros(sz-w_len,1,'single');
%         G_f = zeros(sz-w_len,1,'single');
%         
%         temp2 = g_func_env_tst(g_st,v_avg,r,sz,...
%             b_vl);
%         
%         temp2 = temp2.*10^(-9);
%         
%         %%{
%         for k = 1:length(temp2)
%             temp2(k) = temp2(k) + ((-1).^(round(rand*100)))* temp2(k) *0.10 * rand;
%         end
%         %}
%         %figure(1); plot(temp2(2:end));
%     
%     
%         E_f = temp2*R*W;
%         
%         idx_tmp = single(find(temp2(1:w_len)==max(temp2(1:w_len))));
%         S_idx = idx_tmp(1);
%         tm = [single((S_idx)/100) ...
%             single(linspace(((w_len)/100+0.01),t_end,sz-w_len-1))];
%         
%         E(1) = mean(E_f(1:w_len+1));
%         E(2:(sz-w_len)) = E_f((w_len+2):end);
%         
%         E = log(E);
%         
%         %%{
%         err_sum = g_func_wrpr_mtrx(g_rng,w_len,E,...
%             tm,v_avg,r,sz);
%         
%         [~,g_idx] = mink(err_sum,1);
%         
%         if (g_idx < 100) && (g_idx > 1)
%             
%             g_rng2 = logspace(log10(g_rng(g_idx-1)),...
%                 log10(g_rng(g_idx+1)),100);
%             
%             err_sum = g_func_wrpr_mtrx(g_rng2,w_len,E,...
%                 tm,v_avg,r,sz);
%             
%             [~,g_idx2] = mink(err_sum,1);
%             
%             g_st_w = g_rng2(g_idx2);
%             
%         else
%             
%             g_st_w = g_rng(g_idx);
%         end
%     
%         %}
%         
%         temp3 = g_func_env_tst(g_st_w,v_avg,r,sz,...
%             0);
%         
%         temp1 = temp3.*10^(-9);
%         
%         G_f(1) = mean(temp1(1:w_len+1));
%         G_f(2:(sz-w_len)) = temp1((w_len+2):end);
%         
%         G_f = log(G_f);
%         
%         X = [ones((sz-w_len),1,'single') -tm'];
%         
%         b = X\(E-G_f);
%         
%         E_calc = G_f - b(2).*tm' + b(1);
%         
%         sum_tmp = sum(abs((E(2:end)-E_calc(2:end))))/length(E) * 100;
%         err_tmp = abs(exp(E(1)) - exp(E_calc(1)));
%         error_prom = abs(err_tmp/exp(E(1)))*100;
%         
%         if(b(2) > 0) && (error_prom < 1)
%             
%             b_inv{i} = single(b(2)/(2*pi*freq));
%             g_st_inv{i} = single(g_st_w*v_avg/(2*pi*freq));
%             
%             P_bin12{i} = single(P_d'/sum(P_d));
%             g_orig{i} = single(g_st);
%             b_orig{i} = single(b_vl);
%         else
%             b_inv{i} = single(b_inv{i});
%             g_st_inv{i} = single(g_st_inv{i});
%             
%             P_bin12{i} = single(P_bin12{i});
%             g_orig{i} = single(g_orig{i});
%             b_orig{i} = single(b_orig{i});
%         end
%         %plot(E_calc); hold on; plot(E); hold off;
%         
%         if (rem(i,100)==0)
%         waitbar(i/length(A{1}),f,sprintf('%12.0f',i))
%         end
%         %{
%         g_st
%         g_st_w
%         b_vl
%         b(2)
%         log(W*R)
%         b(1)
%         error_prom
%         %}
%         
%         %{
%                 msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
%                 msft_sgn = sum(E-E_calc);
%                 if (msft_sgn > 0)
%                     E_calc = E_calc + msft;
%                 else
%                     E_calc = E_calc - msft;
%                 end
%                 t_fnl = linspace(-1.0,tm(end),sz-1);
%                 
%                 h1 = figure(j*15); set(gca,'FontName','Helvetica',...
%                     'Fontsize',14,'FontWeight','bold'); hold on;
%                 scatter(tm(1),E(1),'r','LineWidth',2); ...
%                     plot(tm,E,'r','LineWidth',2); %legend(string(r));
%                 hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2);...
%                     plot(tm,E_calc,'b','LineWidth',2); ...
%                     plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
%                     ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}',...
%                     'S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
%                 
%                 figure(j*150); set(gca,'FontName','Helvetica','Fontsize',...
%                     14,'FontWeight','bold'); hold on;...
%                     xlabel('Time (in s)');...
%                     ylabel('Amplitude');
%                 plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.');...
%                     hold off;
%         %}
%     else
%         b_inv{i} = single(b_inv{i});
%         g_st_inv{i} = single(g_st_inv{i});
%         
%         P_bin12{i} = single(P_bin12{i});
%         g_orig{i} = single(g_orig{i});
%         b_orig{i} = single(b_orig{i});
%     end
%     
% end
% 
% P_bin12 = cell2mat(P_bin12);
% g_st_inv = cell2mat(g_st_inv);
% b_inv = cell2mat(b_inv);
% g_orig = cell2mat(g_orig);
% b_orig = cell2mat(b_orig);
% 
% toc
% 
% % delete(f);
% 
% %{
% g_f = lsqlin(double(P_bin12),double((g_st_inv)),[],[],[],[],...
%     ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);
% b_f = lsqlin(double(P_bin12),double((b_inv)),[],[],[],[],...
%     ones(n_grd,1)*1e-5,ones(n_grd,1)*1e-1,[],options);
% 
%  %%%%%%%%%%%%%%%%%%%%%% Or %%%%%%%%%%%%%%%%%%%
% %}
% 
% %{
% 
% g_f = P_bin12\g_st_inv;
% 
% b_f = P_bin12\b_inv;
% 
% idx = find(g_f > 0 & b_f > 0);
% 
% arset = [bincord(idx,:) log10(g_f(idx)) log10(b_f(idx))];
% 
% %}
% 
% %%{
% g_f = sum(bsxfun(@times,g_st_inv,P_bin12))...
%     ./sum(P_bin12);
% b_f = sum(bsxfun(@times,b_inv,P_bin12))...
%     ./sum(P_bin12);
% 
% %arset = [bincord g_f'-10.^G_12{3} b_f'-10.^B_12{3}];
% arset = [bincord g_f' b_f'];
% %}
% 
% save('vars_orig_chkbrd1632.mat','b_orig','g_orig');
% save('vars_inv_chkbrd1632.mat','b_inv','g_st_inv');
% save('prob_chkbrd1632.mat','P_bin12');
% 
% sp = sprintf('%s/bin_chk_31_46_1_16-32Hz.txt',wkdir);
% fid = fopen(sp,'w');
% fprintf(fid,'%f %f %f %f %f\n',arset');
% fclose(fid);
% %}
% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F);
% 
% diffb = abs(b_orig - b_inv)./b_orig * 100;
% diffg = abs(g_orig - g_st_inv)./g_orig * 100;
% figure(1); histogram(diffg);
% figure(2); histogram(diffb);