%% Assigning search directory for events and reading station coordinates and average station velocities 
clc;
close all;
clearvars;


% Directory list for events
dr_lst = [%dir("/home/user/env_EGEL/husn_events/2010*"); dir("/home/user/env_EGEL/husn_events/2011*");...
    %dir("/home/user/env_EGEL/husn_events/2010-09*"); dir("/home/user/env_EGEL/husn_events/2010-10*");
    %dir("/home/user/env_EGEL/husn_events/2010-11*");... %dir("/home/user/env_EGEL/husn_events/2010-12*");
    dir("/home/user/env_EGEL/husn_events/2010*"); dir("/home/user/env_EGEL/husn_events/2011*");...
    dir("/home/user/env_EGEL/husn_events/2012*"); dir("/home/user/env_EGEL/husn_events/2013*"); ...
    dir("/home/user/env_EGEL/husn_events/2014*"); dir("/home/user/env_EGEL/husn_events/2015*"); ...
    dir("/home/user/env_EGEL/husn_events/2016*"); ...
    dir("/home/user/env_EGEL/husn_events/2017*"); ... 
    dir("/home/user/env_EGEL/husn_events/2018*"); ...
    dir("/home/user/env_EGEL/husn_events/2019*"); dir("/home/user/env_EGEL/EGEL_events/200*");...
    dir("/home/user/env_EGEL/Cycnet_events/0*");];

% file for station coordinates fmt= 'Name' 'Latitude' 'Longitude'
fid = fopen('stlist_f.txt','r');
st_db = textscan(fid,'%s %f %f');
fclose(fid);

% file for station velocities fmt= 'Name' 'Avg. Velocity'
fid = fopen('st_vel.txt','r');
velst = textscan(fid,'%s %f');
fclose(fid);

wkdir = '/home/user/env_EGEL/env_gauss2D_chktst';

mkdir(wkdir)

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

b_set = zeros(length(Nx),1);
g_set = zeros(length(Nx),1);

%{

m=1;

for i=1:np_la
    for j=1:np_lo 
        flag = (-1)^(i+j);
        if (flag==1)
            b_set(m) = 0.16;
            g_set(m) = 0.047;
            m = m + 1;
        else
            b_set(m) = 0.262;
            g_set(m) = 0.085;
            m = m + 1;
        end
    end
end
%}


%%{
xgv = 20.4:0.80:30.0;
ygv = 33.4:0.80:39.8;

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
    b_set(In1) = 0.05;
    g_set(In1) = 0.01;
else
    b_set(In1) = 0.005;
    g_set(In1) = 0.001;
end

end
end
%}

chkm_dir = sprintf('%s/bin_chk_orig_bg.txt',wkdir);
fid = fopen(chkm_dir,'w');
fprintf(fid,'%f %f %f %f %f\n',[bincord log10(g_set) log10(b_set)]');
fclose(fid);


%% considering 100 Hz sampling frequency
% Smoothing window size
sm_win12 = 201; sm_win24 = 201; sm_win48 = 201; sm_win816 = 201; sm_win1632 = 201;

% Minimum total length of envelope in samples considering 100 Hz sampling frequency
envl_min = 1900;

% Direct S-wave window length in samples
d_S_len = 1100;


%% Computing for 1-2 Hz band

%abs_pdir = sprintf('%s/ray_prmtrs12.txt',wkdir);
fid_rp = fopen('ray_prmtrs12.txt','r');
A_ry = textscan(fid_rp,"%f %f %f %f %f %f %f %f %f");
fclose(fid_rp);

g_bin12 = zeros(180000,length(Nx),'single');
b_bin12 = zeros(180000,length(Nx),'single');

sum_Pd = zeros(length(Nx),1,'single');

freq = 1.5;

ri = single(1);


%w_len = single(d_S_len);

for i=1:13%length(dr_lst)
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_1-2.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        n_fls = length(aw);
        %n_fls = 1;
        if (n_fls > 0)
            srch_str = sprintf('%s/*.hyp',chek);
            [hyp_dir,hyp_info] = grep('-s',"GEOGRAPHIC",srch_str);
            evnt_info = strsplit(hyp_info.match{hyp_info.nfiles},' ');
            ev_lat = str2double(evnt_info{10});
            ev_lon = str2double(evnt_info{12});
            ev_dep = str2double(evnt_info{14});
            
            Earset_cl = cell(1,n_fls);
            r = zeros(n_fls,1);
            v_avg = zeros(n_fls,1);
            S_idx = zeros(n_fls,1);
            sz = zeros(n_fls,1);
            t_end = zeros(n_fls,1);
            g_st = zeros(n_fls,1);
            b_vl = zeros(n_fls,1);
            tm_cl = cell(1,n_fls);
            stnm = cell(1,n_fls);
            
            st_lat = zeros(n_fls,1);
            st_lon = zeros(n_fls,1);
            
            [ev_x,ev_y,ev_z] = geodetic2ecef(wgs84,ev_lat,ev_lon,ev_dep);
            
            W = rand*10000;
            R = rand(n_fls,1)*10;
            R = R/geomean(R);
            
            for j = 1:n_fls
            
                q1 = split(aw(j).name,'_');
                stnm{j} = q1(3);
                vel_fl = find(strcmp(velst{1},stnm{j}));
                if (vel_fl)
                    v_avg(j) = single(velst{2}(vel_fl));
                else
                    v_avg(j) = single(3300);
                end
                idx_st = find(strcmp(st_db{1},stnm{j}));
                st_lat(j) = st_db{2}(idx_st);
                st_lon(j) = st_db{3}(idx_st);
                [st_x,st_y,st_z] = geodetic2ecef(wgs84,st_lat(j),st_lon(j),0);
                r(j) = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
                sz(j) = aw(j).bytes;
                sz(j) = round(sz(j)/9);
                t_end(j) = (sz(j)-100)/100;
                %r(j)
                %v_avg(j)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                s_mj = (v_avg(j)/1000 * t_end(j) + r(j))/2;
                s_mn = sqrt(s_mj^2 - (r(j)/2)^2);
                Cx = (st_x + ev_x)/2;
                Cy = (st_y + ev_y)/2;
                slp_ry = (st_y-ev_y)/(st_x-ev_x);
                dst_ry = abs((Ny-Cy) - slp_ry*(Nx-Cx))./sqrt(1+slp_ry^2);
                slp_rypp = -1/slp_ry;
                dst_rypp = abs((Ny-Cy) - slp_rypp*(Nx-Cx))./sqrt(1+slp_rypp^2);
                
                P_sgma = [s_mj^2/9 s_mn^2/9];
                P_d = mvnpdf([dst_rypp dst_ry],[],P_sgma);
                P_d = normc(P_d);
                sum_Pd = sum_Pd + P_d;
                
                dst = (dst_rypp).^2./(s_mj)^2 + (dst_ry).^2./(s_mn)^2;
                
                idx = find(dst < 1);
                
                tempg = g_set(idx) .* P_d(idx);
                tempb = b_set(idx) .* P_d(idx);
                
                s_g = length(tempg);
                
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
                g_st(j) = mean_g;
                b_vl(j) = mean_b;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                w_len=1100;
                
                r(j) = r(j)*1000;
                
                temp2 = g_func_env_tst(g_st(j),v_avg(j),r(j),sz(j),w_len-100,t_end(j),b_vl(j));
                
                %%{
                for k = 1:length(temp2)
                    temp2(k) = temp2(k) + ((-1).^(round(rand*100)))* temp2(k) *0.1 * rand;
                end
                %}
                
                E_f = log(temp2*R(j)*W);
                 
                idx_tmp = single(find(temp2(1:w_len)==max(temp2(1:w_len))));
                S_idx(j) = idx_tmp(1);
                tme = [single((S_idx(j))/100) single(linspace(((w_len-100)/100+0.01),t_end(j),sz(j)-1-w_len))];
                tm_cl{j} = tme;
                
                Earset_cl{j}(1) = mean(E_f(1:w_len));
                Earset_cl{j}(2:(sz(j)-w_len)) = E_f((w_len+1):end-1);
                %length(tme)
                %length(Earset_cl{j})
            end
            
            Earset = cell2mat(Earset_cl);
            tm = cell2mat(tm_cl);
            
            
            opts = optimoptions('lsqnonlin');
            %opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            opts.Display = 'off';
            g_sti = ones(1,n_fls) * 1e-5;
            sum_tmp = 0;
            error_prom = 0;
            err_tmp = 0;
            
            R2 = rand(n_fls,1)*10;
            R2 = R2/geomean(R2);
            
            fun = @(g_sti)g_func_wrapper_weg(g_sti,w_len,Earset_cl,tm_cl,v_avg,r,sz,n_fls);
            g_st_inv = lsqnonlin(fun,g_sti,zeros(n_fls,1),ones(n_fls,1),opts);
            
            E = single(Earset);
            G_f_set = cell(1,n_fls);
            E_cl = cell(1,n_fls);
            Cnst = cell(1,n_fls);
            b_inv = zeros(n_fls,2);
            for j=1:n_fls
                v = v_avg(j);
                r_d = r(j);
                sz_d = sz(j);
                
                t_end2 = (sz_d-100)/100;
                
                temp1 = g_func_env_tst(g_st_inv(j),v,r_d,sz_d,(w_len-100),t_end2,0);
                temp3 = log(temp1);
                G_f_set{j}(1) = mean(temp3(1:w_len));
                G_f_set{j}(2:(sz_d-w_len)) = temp3((w_len+1):end-1);
                Cnst{j} = ones(1,(sz_d-w_len),'single') * R2(j);
                X = [ones((sz_d-w_len),1,'single') -tm_cl{j}'];
                lb = [0 1e-6];
                %ub = [10 0.5];
                options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
                options.Display = 'off';
                options = optimoptions(options,'SubproblemAlgorithm','factorization');
                %options = optimoptions(options,'FunctionTolerance',1.000000e-10,'OptimalityTolerance',1e-34);
                [b,resnorm,residual,exitflag,output] = ...
                    lsqlin(double(X),double((Earset_cl{j}-G_f_set{j})'),[],[],[],[],lb,[],[],options);
                b_inv(j,:) = b;
                E_cl{j} = G_f_set{j}-b(2).*tm_cl{j} + b(1);
                
                sum_tmp = sum_tmp + sum(abs((Earset_cl{j}-E_cl{j}))./E_cl{j});
                err_tmp = err_tmp + abs((E_cl{j}(1) - Earset_cl{j}(1)));
                error_prom = error_prom + abs(err_tmp/Earset_cl{j}(1));
            end
            %{
            G_f = cell2mat(G_f_set);
            Cnst_ar = cell2mat(Cnst);
            
            X = [ones(length(G_f),1,'single') -tm'];

            lb = [0 0];
            %ub = [10 0.5];
            options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
            options.Display = 'off';
            %options = optimoptions(options,'SubproblemAlgorithm','factorization');
            %options = optimoptions(options,'FunctionTolerance',1.000000e-10,'OptimalityTolerance',1e-34);
            [b_inv,resnorm,residual,exitflag,output] = ...
                lsqlin(double(X),double((E-G_f)'),[],[],[],[],lb,[],[],options);
            %b = X\(E-G_f)';
            E_calc = G_f-b_inv(2).*tm + b_inv(1);% + b_inv(2);
            %}
            E_calc = cell2mat(E_cl);
            %plot(E_calc); hold on; plot(Earset); hold off;
            g_st
            g_st_inv
            b_vl
            b_inv
            %R2*b_inv(2)
            %R
            W
            log(W)
            mean(b_inv(:,1))
            %if(~isnan(b_inv(:,2))) && (b_inv(:,2) > 0)
                
                if (err_tmp <= 1.18) && (sum_tmp <= 0.50)
                    %disp('success');
                    for j=1:n_fls
                    g_bin12(ri,idx) = P_d(idx) * g_st_inv(j)/sum(P_d);
                    
                    b_bin12(ri,idx) = P_d(idx) * b_inv(j,2)/sum(P_d);
                    ri = ri + 1;
                    end
                end
                %{
        msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
        msft_sgn = sum(E-E_calc);
        if (msft_sgn > 0)
            E_calc = E_calc + msft;
        else
            E_calc = E_calc - msft;
        end
        t_fnl = linspace(-1.0,tm(end),sz-1);
        
        h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
        scatter(tm(1),E(1),'r','LineWidth',2); plot(tm,E,'r','LineWidth',2); %legend(string(r));
        hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2); plot(tm,E_calc,'b','LineWidth',2); ...
            plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
            ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
        
        figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
            xlabel('Time (in s)');...
            ylabel('Amplitude');
        plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.'); hold off;
                %}
            %end
        end
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

arset = [bincord log10(g_f) log10(b_f)];

sp = sprintf('%s/bin_chk_31_46_1_1-2Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f %f\n',arset');
fclose(fid);
%{

%% Computing for 2-4 Hz band


freq = 3;

count = zeros(length(Nx),1);
g_avg = zeros(length(Nx),1);
b_avg = zeros(length(Nx),1);

sum_Pd = 0;

%w_len = single(d_S_len);
flag = 0;

for i=1:length(dr_lst)
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_2-4.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        n_fls = length(aw);
        if (n_fls > 0)
            srch_str = sprintf('%s/*.hyp',chek);
            [hyp_dir,hyp_info] = grep('-s',"GEOGRAPHIC",srch_str);
            evnt_info = strsplit(hyp_info.match{hyp_info.nfiles},' ');
            ev_lat = str2double(evnt_info{10});
            ev_lon = str2double(evnt_info{12});
            ev_dep = str2double(evnt_info{14});
            
            b_inv = zeros(2,n_fls,'single');
            b = zeros(2,n_fls,'single');
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            opts.Display = 'off';
            g_sti = 1e-4;
            g_st_inv = ones(1,n_fls) * 1e-4;
            g_st = ones(1,n_fls) * 1e-4;
            sum_tmp = 0;
            error_prom = 0; 
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,aw(j).name);
                q1 = split(aw(j).name,'_');
                stnm = q1(3);
                
                vel_fl = find(strcmp(velst{1},stnm));
                if (vel_fl)
                    v_avg = single(velst{2}(vel_fl));
                else
                    v_avg = single(3300);
                end
                v_avg = v_avg/1000;
                idx_st = find(strcmp(st_db{1},stnm));
                [ev_x,ev_y,ev_z] = geodetic2ecef(wgs84,ev_lat,ev_lon,ev_dep);
                [st_x,st_y,st_z] = geodetic2ecef(wgs84,st_db{2}(idx_st),st_db{3}(idx_st),0);
                r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
                sz = aw(j).bytes;
                sz = round(sz/9);
                t_end = sz/100;
                s_mj = (v_avg * t_end + r)/2;
                s_mn = sqrt(s_mj^2 - (r/2)^2);
                Cx = (st_x + ev_x)/2;
                Cy = (st_y + ev_y)/2;
                Cz = (st_z + ev_z)/2;
                dst = ((Nx - Cx).^2)./(s_mj^2) + ((Ny - Cy).^2)./(s_mn^2);
                P_sgma = [s_mj/3 s_mn/3];
                P_mu = [Cx Cy];
                P_d = mvnpdf([Nx Ny],P_mu,P_sgma);
                sum_Pd = P_d + sum_Pd;
                pos_idx = find(P_d > 0);
                g_st(j) = sum(g_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                b(j) = sum(b_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                w_len=1100;
                if (sz > 1900)
                    temp2 = g_func_env_tst(g_st(j),v_avg,r,sz,w_len-100,t_end,b(j));
                    E_f = temp2;
                    %E_f = log(temp1);
                    %E(1)
                    %{
                    for k = 1:length(E_f)
                        E_f(k) = E_f(k) + (-1).^(round(rand*100))* E_f(k) * 0.1 * rand;
                    end
                    %}
                    %E(1)
                    S_idx = single(find(E_f(1:w_len)==max(E_f(1:w_len))));
                    E = ones(1,(sz-w_len),'single');
                    
                    E(1) = sum(E_f(1:w_len));
                    E(2:(sz-w_len)) = E_f((w_len+1):end-1);
                    E = log(E);
                    g_sti = 1e-4;
                    tm = zeros(sz-w_len,1,'single');
                    tm(1) = S_idx(1)/100;
                    %tm(1) = (median(find(abs(log(E_arr(1:end))-E(1)) < 0.01))-101)/100;
                    tm(2:sz-w_len) = linspace(((w_len-101)/100+0.01),t_end,sz-w_len-1);
                    fun = @(g_sti)g_func_wrapper2(g_sti,w_len,E,tm,v_avg,r,t_end);
                    g_st_inv(j) = lsqnonlin(fun,g_sti,[],[],opts);
                    temp1 = g_func_env(g_st(j),v_avg,r,(sz-w_len),(w_len-100),t_end);
                    G_f = log(temp1);
                    
                    X = [ones(length(G_f),1) -tm];
                    b_inv(:,j) = double(X)\(double(E-G_f))';
                    %g_st_inv(j)
                    %g_st(j)
                    if(~isnan(b_inv(2,j)))
                        E_calc = G_f+b_inv(1,j)-b_inv(2,j).*tm';
                        
                        sum_tmp = sum_tmp + sum(abs((E-E_calc)/E));
                        err_tmp = abs((E_calc(1) - E(1)));
                        error_prom = error_prom + abs(err_tmp/E(1));
                        
                        if (err_tmp < 1)
                            %disp('success');
                            g_inv_set = (g_st(j) * P_d);
                            b_inv_set = (b_inv(2,j) * P_d);
                            count(pos_idx) = count(pos_idx) + 1;
                            g_avg(pos_idx) = (g_avg(pos_idx) + g_inv_set(pos_idx));
                            b_avg(pos_idx) = (b_avg(pos_idx) + b_inv_set(pos_idx));
                        end
                        %{
                        msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                        msft_sgn = sum(E-E_calc);
                        if (msft_sgn > 0)
                            E_calc = E_calc + msft;
                        else
                            E_calc = E_calc - msft;
                        end
                        t_fnl = linspace(-1.0,tm(end),sz-1);
                        
                        h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                        scatter(tm(1),E(1),'r','LineWidth',2); plot(tm,E,'r','LineWidth',2); %legend(string(r));
                        hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2); plot(tm,E_calc,'b','LineWidth',2); ...
                            plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                            ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                        
                        figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                            xlabel('Time (in s)');...
                            ylabel('Amplitude');
                        plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.'); hold off;
                        %}
                    end
                end
            end
        end
        if (flag==1)
            break;
        end
    end
end

arset = [bincord log10(g_avg./sum_Pd) log10(b_avg./sum_Pd) count];

sp = sprintf('%s/bin_chk_31_46_1_2-4Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f %f\n',arset');
fclose(fid);

%% Computing for 4-8 Hz

freq = 6;

count = zeros(length(Nx),1);
g_avg = zeros(length(Nx),1);
b_avg = zeros(length(Nx),1);

sum_Pd = 0;

%w_len = single(d_S_len);
flag = 0;

for i=1:length(dr_lst)
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_4-8.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        n_fls = length(aw);
        if (n_fls > 0)
            srch_str = sprintf('%s/*.hyp',chek);
            [hyp_dir,hyp_info] = grep('-s',"GEOGRAPHIC",srch_str);
            evnt_info = strsplit(hyp_info.match{hyp_info.nfiles},' ');
            ev_lat = str2double(evnt_info{10});
            ev_lon = str2double(evnt_info{12});
            ev_dep = str2double(evnt_info{14});
            
            b_inv = zeros(2,n_fls,'single');
            b = zeros(2,n_fls,'single');
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            opts.Display = 'off';
            g_sti = 1e-4;
            g_st_inv = ones(1,n_fls) * 1e-4;
            g_st = ones(1,n_fls) * 1e-4;
            sum_tmp = 0;
            error_prom = 0; 
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,aw(j).name);
                q1 = split(aw(j).name,'_');
                stnm = q1(3);
                
                vel_fl = find(strcmp(velst{1},stnm));
                if (vel_fl)
                    v_avg = single(velst{2}(vel_fl));
                else
                    v_avg = single(3300);
                end
                v_avg = v_avg/1000;
                idx_st = find(strcmp(st_db{1},stnm));
                [ev_x,ev_y,ev_z] = geodetic2ecef(wgs84,ev_lat,ev_lon,ev_dep);
                [st_x,st_y,st_z] = geodetic2ecef(wgs84,st_db{2}(idx_st),st_db{3}(idx_st),0);
                r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
                sz = aw(j).bytes;
                sz = round(sz/9);
                t_end = sz/100;
                s_mj = (v_avg * t_end + r)/2;
                s_mn = sqrt(s_mj^2 - (r/2)^2);
                Cx = (st_x + ev_x)/2;
                Cy = (st_y + ev_y)/2;
                Cz = (st_z + ev_z)/2;
                dst = ((Nx - Cx).^2)./(s_mj^2) + ((Ny - Cy).^2)./(s_mn^2);
                P_sgma = [s_mj/3 s_mn/3];
                P_mu = [Cx Cy];
                P_d = mvnpdf([Nx Ny],P_mu,P_sgma);
                sum_Pd = P_d + sum_Pd;
                pos_idx = find(P_d > 0);
                g_st(j) = sum(g_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                b(j) = sum(b_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                w_len=1100;
                if (sz > 1900)
                    temp2 = g_func_env_tst(g_st(j),v_avg,r,sz,w_len-100,t_end,b(j));
                    E_f = temp2;
                    %E_f = log(temp1);
                    %E(1)
                    %{
                    for k = 1:length(E_f)
                        E_f(k) = E_f(k) + (-1).^(round(rand*100))* E_f(k) * 0.1 * rand;
                    end
                    %}
                    %E(1)
                    S_idx = single(find(E_f(1:w_len)==max(E_f(1:w_len))));
                    E = ones(1,(sz-w_len),'single');
                    
                    E(1) = sum(E_f(1:w_len));
                    E(2:(sz-w_len)) = E_f((w_len+1):end-1);
                    E = log(E);
                    g_sti = 1e-4;
                    tm = zeros(sz-w_len,1,'single');
                    tm(1) = S_idx(1)/100;
                    %tm(1) = (median(find(abs(log(E_arr(1:end))-E(1)) < 0.01))-101)/100;
                    tm(2:sz-w_len) = linspace(((w_len-101)/100+0.01),t_end,sz-w_len-1);
                    fun = @(g_sti)g_func_wrapper2(g_sti,w_len,E,tm,v_avg,r,t_end);
                    g_st_inv(j) = lsqnonlin(fun,g_sti,[],[],opts);
                    temp1 = g_func_env(g_st(j),v_avg,r,(sz-w_len),(w_len-100),t_end);
                    G_f = log(temp1);
                    
                    X = [ones(length(G_f),1) -tm];
                    b_inv(:,j) = double(X)\(double(E-G_f))';
                    %g_st_inv(j)
                    %g_st(j)
                    if(~isnan(b_inv(2,j)))
                        E_calc = G_f+b_inv(1,j)-b_inv(2,j).*tm';
                        
                        sum_tmp = sum_tmp + sum(abs((E-E_calc)/E));
                        err_tmp = abs((E_calc(1) - E(1)));
                        error_prom = error_prom + abs(err_tmp/E(1));
                        
                        if (err_tmp < 1)
                            %disp('success');
                            g_inv_set = (g_st(j) * P_d);
                            b_inv_set = (b_inv(2,j) * P_d);
                            count(pos_idx) = count(pos_idx) + 1;
                            g_avg(pos_idx) = (g_avg(pos_idx) + g_inv_set(pos_idx));
                            b_avg(pos_idx) = (b_avg(pos_idx) + b_inv_set(pos_idx));
                        end
                        %{
                        msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                        msft_sgn = sum(E-E_calc);
                        if (msft_sgn > 0)
                            E_calc = E_calc + msft;
                        else
                            E_calc = E_calc - msft;
                        end
                        t_fnl = linspace(-1.0,tm(end),sz-1);
                        
                        h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                        scatter(tm(1),E(1),'r','LineWidth',2); plot(tm,E,'r','LineWidth',2); %legend(string(r));
                        hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2); plot(tm,E_calc,'b','LineWidth',2); ...
                            plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                            ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                        
                        figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                            xlabel('Time (in s)');...
                            ylabel('Amplitude');
                        plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.'); hold off;
                        %}
                    end
                end
            end
        end
        if (flag==1)
            break;
        end
    end
end

arset = [bincord log10(g_avg./sum_Pd) log10(b_avg./sum_Pd) count];

sp = sprintf('%s/bin_chk_31_46_1_4-8Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f %f\n',arset');
fclose(fid);
%% Computing for 8-16 Hz

freq = 12;

count = zeros(length(Nx),1);
g_avg = zeros(length(Nx),1);
b_avg = zeros(length(Nx),1);

sum_Pd = 0;

%w_len = single(d_S_len);
flag = 0;

for i=1:length(dr_lst)
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_8-16.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        n_fls = length(aw);
        if (n_fls > 0)
            srch_str = sprintf('%s/*.hyp',chek);
            [hyp_dir,hyp_info] = grep('-s',"GEOGRAPHIC",srch_str);
            evnt_info = strsplit(hyp_info.match{hyp_info.nfiles},' ');
            ev_lat = str2double(evnt_info{10});
            ev_lon = str2double(evnt_info{12});
            ev_dep = str2double(evnt_info{14});
            
            b_inv = zeros(2,n_fls,'single');
            b = zeros(2,n_fls,'single');
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            opts.Display = 'off';
            g_sti = 1e-4;
            g_st_inv = ones(1,n_fls) * 1e-4;
            g_st = ones(1,n_fls) * 1e-4;
            sum_tmp = 0;
            error_prom = 0; 
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,aw(j).name);
                q1 = split(aw(j).name,'_');
                stnm = q1(3);
                
                vel_fl = find(strcmp(velst{1},stnm));
                if (vel_fl)
                    v_avg = single(velst{2}(vel_fl));
                else
                    v_avg = single(3300);
                end
                v_avg = v_avg/1000;
                idx_st = find(strcmp(st_db{1},stnm));
                [ev_x,ev_y,ev_z] = geodetic2ecef(wgs84,ev_lat,ev_lon,ev_dep);
                [st_x,st_y,st_z] = geodetic2ecef(wgs84,st_db{2}(idx_st),st_db{3}(idx_st),0);
                r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
                sz = aw(j).bytes;
                sz = round(sz/9);
                t_end = sz/100;
                s_mj = (v_avg * t_end + r)/2;
                s_mn = sqrt(s_mj^2 - (r/2)^2);
                Cx = (st_x + ev_x)/2;
                Cy = (st_y + ev_y)/2;
                Cz = (st_z + ev_z)/2;
                dst = ((Nx - Cx).^2)./(s_mj^2) + ((Ny - Cy).^2)./(s_mn^2);
                P_sgma = [s_mj/3 s_mn/3];
                P_mu = [Cx Cy];
                P_d = mvnpdf([Nx Ny],P_mu,P_sgma);
                sum_Pd = P_d + sum_Pd;
                pos_idx = find(P_d > 0);
                g_st(j) = sum(g_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                b(j) = sum(b_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                w_len=1100;
                if (sz > 1900)
                    temp2 = g_func_env_tst(g_st(j),v_avg,r,sz,w_len-100,t_end,b(j));
                    E_f = temp2;
                    %E_f = log(temp1);
                    %E(1)
                    %{
                    for k = 1:length(E_f)
                        E_f(k) = E_f(k) + (-1).^(round(rand*100))* E_f(k) * 0.1 * rand;
                    end
                    %}
                    %E(1)
                    S_idx = single(find(E_f(1:w_len)==max(E_f(1:w_len))));
                    E = ones(1,(sz-w_len),'single');
                    
                    E(1) = sum(E_f(1:w_len));
                    E(2:(sz-w_len)) = E_f((w_len+1):end-1);
                    E = log(E);
                    g_sti = 1e-4;
                    tm = zeros(sz-w_len,1,'single');
                    tm(1) = S_idx(1)/100;
                    %tm(1) = (median(find(abs(log(E_arr(1:end))-E(1)) < 0.01))-101)/100;
                    tm(2:sz-w_len) = linspace(((w_len-101)/100+0.01),t_end,sz-w_len-1);
                    fun = @(g_sti)g_func_wrapper2(g_sti,w_len,E,tm,v_avg,r,t_end);
                    g_st_inv(j) = lsqnonlin(fun,g_sti,[],[],opts);
                    temp1 = g_func_env(g_st(j),v_avg,r,(sz-w_len),(w_len-100),t_end);
                    G_f = log(temp1);
                    
                    X = [ones(length(G_f),1) -tm];
                    b_inv(:,j) = double(X)\(double(E-G_f))';
                    %g_st_inv(j)
                    %g_st(j)
                    if(~isnan(b_inv(2,j)))
                        E_calc = G_f+b_inv(1,j)-b_inv(2,j).*tm';
                        
                        sum_tmp = sum_tmp + sum(abs((E-E_calc)/E));
                        err_tmp = abs((E_calc(1) - E(1)));
                        error_prom = error_prom + abs(err_tmp/E(1));
                        
                        if (err_tmp < 1)
                            %disp('success');
                            g_inv_set = (g_st(j) * P_d);
                            b_inv_set = (b_inv(2,j) * P_d);
                            count(pos_idx) = count(pos_idx) + 1;
                            g_avg(pos_idx) = (g_avg(pos_idx) + g_inv_set(pos_idx));
                            b_avg(pos_idx) = (b_avg(pos_idx) + b_inv_set(pos_idx));
                        end
                        %{
                        msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                        msft_sgn = sum(E-E_calc);
                        if (msft_sgn > 0)
                            E_calc = E_calc + msft;
                        else
                            E_calc = E_calc - msft;
                        end
                        t_fnl = linspace(-1.0,tm(end),sz-1);
                        
                        h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                        scatter(tm(1),E(1),'r','LineWidth',2); plot(tm,E,'r','LineWidth',2); %legend(string(r));
                        hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2); plot(tm,E_calc,'b','LineWidth',2); ...
                            plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                            ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                        
                        figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                            xlabel('Time (in s)');...
                            ylabel('Amplitude');
                        plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.'); hold off;
                        %}
                    end
                end
            end
        end
        if (flag==1)
            break;
        end
    end
end

arset = [bincord log10(g_avg./sum_Pd) log10(b_avg./sum_Pd) count];

sp = sprintf('%s/bin_chk_31_46_1_8-16Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f %f\n',arset');
fclose(fid);

%% Computing for 16-32 Hz band

freq = 24;

count = zeros(length(Nx),1);
g_avg = zeros(length(Nx),1);
b_avg = zeros(length(Nx),1);

sum_Pd = 0;

%w_len = single(d_S_len);
flag = 0;

for i=1:length(dr_lst)
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_16-32.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        n_fls = length(aw);
        if (n_fls > 0)
            srch_str = sprintf('%s/*.hyp',chek);
            [hyp_dir,hyp_info] = grep('-s',"GEOGRAPHIC",srch_str);
            evnt_info = strsplit(hyp_info.match{hyp_info.nfiles},' ');
            ev_lat = str2double(evnt_info{10});
            ev_lon = str2double(evnt_info{12});
            ev_dep = str2double(evnt_info{14});
            
            b_inv = zeros(2,n_fls,'single');
            b = zeros(2,n_fls,'single');
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            opts.Display = 'off';
            g_sti = 1e-4;
            g_st_inv = ones(1,n_fls) * 1e-4;
            g_st = ones(1,n_fls) * 1e-4;
            sum_tmp = 0;
            error_prom = 0; 
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,aw(j).name);
                q1 = split(aw(j).name,'_');
                stnm = q1(3);
                
                vel_fl = find(strcmp(velst{1},stnm));
                if (vel_fl)
                    v_avg = single(velst{2}(vel_fl));
                else
                    v_avg = single(3300);
                end
                v_avg = v_avg/1000;
                idx_st = find(strcmp(st_db{1},stnm));
                [ev_x,ev_y,ev_z] = geodetic2ecef(wgs84,ev_lat,ev_lon,ev_dep);
                [st_x,st_y,st_z] = geodetic2ecef(wgs84,st_db{2}(idx_st),st_db{3}(idx_st),0);
                r = sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2);
                sz = aw(j).bytes;
                sz = round(sz/9);
                t_end = sz/100;
                s_mj = (v_avg * t_end + r)/2;
                s_mn = sqrt(s_mj^2 - (r/2)^2);
                Cx = (st_x + ev_x)/2;
                Cy = (st_y + ev_y)/2;
                Cz = (st_z + ev_z)/2;
                dst = ((Nx - Cx).^2)./(s_mj^2) + ((Ny - Cy).^2)./(s_mn^2);
                P_sgma = [s_mj/3 s_mn/3];
                P_mu = [Cx Cy];
                P_d = mvnpdf([Nx Ny],P_mu,P_sgma);
                sum_Pd = P_d + sum_Pd;
                pos_idx = find(P_d > 0);
                g_st(j) = sum(g_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                b(j) = sum(b_set(pos_idx) .* P_d(pos_idx))./sum(P_d(pos_idx));
                w_len=1100;
                if (sz > 1900)
                    temp2 = g_func_env_tst(g_st(j),v_avg,r,sz,w_len-100,t_end,b(j));
                    E_f = temp2;
                    %E_f = log(temp1);
                    %E(1)
                    %{
                    for k = 1:length(E_f)
                        E_f(k) = E_f(k) + (-1).^(round(rand*100))* E_f(k) * 0.1 * rand;
                    end
                    %}
                    %E(1)
                    S_idx = single(find(E_f(1:w_len)==max(E_f(1:w_len))));
                    E = ones(1,(sz-w_len),'single');
                    
                    E(1) = sum(E_f(1:w_len));
                    E(2:(sz-w_len)) = E_f((w_len+1):end-1);
                    E = log(E);
                    g_sti = 1e-4;
                    tm = zeros(sz-w_len,1,'single');
                    tm(1) = S_idx(1)/100;
                    %tm(1) = (median(find(abs(log(E_arr(1:end))-E(1)) < 0.01))-101)/100;
                    tm(2:sz-w_len) = linspace(((w_len-101)/100+0.01),t_end,sz-w_len-1);
                    fun = @(g_sti)g_func_wrapper2(g_sti,w_len,E,tm,v_avg,r,t_end);
                    g_st_inv(j) = lsqnonlin(fun,g_sti,[],[],opts);
                    temp1 = g_func_env(g_st(j),v_avg,r,(sz-w_len),(w_len-100),t_end);
                    G_f = log(temp1);
                    
                    X = [ones(length(G_f),1) -tm];
                    b_inv(:,j) = double(X)\(double(E-G_f))';
                    %g_st_inv(j)
                    %g_st(j)
                    if(~isnan(b_inv(2,j)))
                        E_calc = G_f+b_inv(1,j)-b_inv(2,j).*tm';
                        
                        sum_tmp = sum_tmp + sum(abs((E-E_calc)/E));
                        err_tmp = abs((E_calc(1) - E(1)));
                        error_prom = error_prom + abs(err_tmp/E(1));
                        
                        if (err_tmp < 1)
                            %disp('success');
                            g_inv_set = (g_st(j) * P_d);
                            b_inv_set = (b_inv(2,j) * P_d);
                            count(pos_idx) = count(pos_idx) + 1;
                            g_avg(pos_idx) = (g_avg(pos_idx) + g_inv_set(pos_idx));
                            b_avg(pos_idx) = (b_avg(pos_idx) + b_inv_set(pos_idx));
                        end
                        %{
                        msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                        msft_sgn = sum(E-E_calc);
                        if (msft_sgn > 0)
                            E_calc = E_calc + msft;
                        else
                            E_calc = E_calc - msft;
                        end
                        t_fnl = linspace(-1.0,tm(end),sz-1);
                        
                        h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                        scatter(tm(1),E(1),'r','LineWidth',2); plot(tm,E,'r','LineWidth',2); %legend(string(r));
                        hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2); plot(tm,E_calc,'b','LineWidth',2); ...
                            plot(t_fnl,log(E),'g'); xlabel('Time (in s)');...
                            ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                        
                        figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                            xlabel('Time (in s)');...
                            ylabel('Amplitude');
                        plot(t_fnl,E,'LineWidth',2); legend('Smoothed S-wave env.'); hold off;
                        %}
                    end
                end
            end
        end
        if (flag==1)
            break;
        end
    end
end

arset = [bincord log10(g_avg./sum_Pd) log10(b_avg./sum_Pd) count];

sp = sprintf('%s/bin_chk_31_46_1_16-32Hz.txt',wkdir);
fid = fopen(sp,'w');
fprintf(fid,'%f %f %f %f %f %f\n',arset');
fclose(fid);
%}