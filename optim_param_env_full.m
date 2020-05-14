%% Assigning search directory for events and reading station coordinates
%% and average station velocities 
clc;
close all;
clearvars;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

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
    dir("/home/user/env_EGEL/Cycnet_events/0*"); ...
    dir("/home/user/env_EGEL/kef_events/2*"); dir("/home/user/env_EGEL/cornoth_events/2*");];

% file for station coordinates fmt= 'Name' 'Latitude' 'Longitude'
fid = fopen('stlist_f.txt','r');
st_db = textscan(fid,'%s %f %f');
fclose(fid);

% file for station velocities fmt= 'Name' 'Avg. Velocity'
fid = fopen('st_vel.txt','r');
velst = textscan(fid,'%s %f');
fclose(fid);

wkdir = '/home/user/env_EGEL/env_gauss_6,3,1.5,.75_noavg';

mkdir(wkdir)

wgs84 = wgs84Ellipsoid('kilometers');

%% considering 100 Hz sampling frequency
% Smoothing window size
sm_win12 = 601; sm_win24 = 301; sm_win48 = 151; sm_win816 = 75;
sm_win1632 = 32;

% Minimum total length of envelope in samples considering 100 Hz
% sampling frequency
envl_min = single(1900);

% Direct S-wave window length in samples
d_S_len = single(100);


%% Computing for all frequency bands

abs_pdir = sprintf('%s/ray_prmtrs12.txt',wkdir);
fid_rp12 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_W_1-2.txt',wkdir);
fidW12 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_R_1-2.txt',wkdir);
fidR12 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_err_1-2.txt',wkdir);
fid_er12 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/ray_prmtrs24.txt',wkdir);
fid_rp24 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_W_2-4.txt',wkdir);
fidW24 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_R_2-4.txt',wkdir);
fidR24 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_err_2-4.txt',wkdir);
fid_er24 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/ray_prmtrs48.txt',wkdir);
fid_rp48 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_W_4-8.txt',wkdir);
fidW48 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_R_4-8.txt',wkdir);
fidR48 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_err_4-8.txt',wkdir);
fid_er48 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/ray_prmtrs816.txt',wkdir);
fid_rp816 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_W_8-16.txt',wkdir);
fidW816 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_R_8-16.txt',wkdir);
fidR816 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_err_8-16.txt',wkdir);
fid_er816 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/ray_prmtrs1632.txt',wkdir);
fid_rp1632 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_W_16-32.txt',wkdir);
fidW1632 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_R_16-32.txt',wkdir);
fidR1632 = fopen(abs_pdir,'w');

abs_pdir = sprintf('%s/env_err_16-32.txt',wkdir);
fid_er1632 = fopen(abs_pdir,'w');


f = waitbar(0,'1','Name','Iterating through files for all freq',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

opts = optimoptions('lsqnonlin');
opts.Algorithm = 'levenberg-marquardt';
%opts.MaxIterations = 5000;
opts.Display = 'off';
%opts.StepTolerance=1e-20;


options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
options.Display = 'off';
%options = optimoptions(options,'SubproblemAlgorithm','factorization');

g_rng = single(logspace(-7,1,100));

tic

%%% set i=26212 & j = 60+ for example fit and envelope

for i=26212:26212%length(dr_lst)
    %i
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    %dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_*.bin');
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
            
            [ev_x,ev_y,ev_z] = geodetic2enu(ev_lat,ev_lon,...
                                        ev_dep,34,18,0,wgs84);

            lgW_tmp12 = zeros(n_fls,2,'single');
            lgW_tmp24 = zeros(n_fls,2,'single');
            lgW_tmp48 = zeros(n_fls,2,'single');
            lgW_tmp816 = zeros(n_fls,2,'single');
            lgW_tmp1632 = zeros(n_fls,2,'single');
            
            g_st_tmp = zeros(n_fls,1,'single');
            sum_tmp = zeros(n_fls,1,'single');
            err_tmp = zeros(n_fls,1,'single');
            error_prom = zeros(n_fls,1,'single');
            
            sz = zeros(n_fls,1,'single');
            
            st_lat = zeros(n_fls,1,'single');
            st_lon = zeros(n_fls,1,'single');
            idx_st = zeros(1,n_fls,'single');
            
            b = zeros(2,n_fls,'single');
            g_sti = 1e-4;
            for j = 60:62%n_fls
                fw = split(aw(j).name,'_');
                if (strcmp(fw{5},'1-2.bin'))
                    freq = 1.5;
                    fid_rp = fid_rp12;
                    fid_er = fid_er12;
                    sm_win = sm_win12;
                elseif (strcmp(fw{5},'2-4.bin'))
                    freq = 3;
                    fid_rp = fid_rp24;
                    fid_er = fid_er24;
                    sm_win = sm_win24;
                elseif (strcmp(fw{5},'4-8.bin'))
                    freq = 6;
                    fid_rp = fid_rp48;
                    fid_er = fid_er48;
                    sm_win = sm_win48;
                elseif (strcmp(fw{5},'8-16.bin'))
                    freq = 12;
                    fid_rp = fid_rp816;
                    fid_er = fid_er816;
                    sm_win = sm_win816;
                else
                    freq = 24;
                    fid_rp = fid_rp1632;
                    fid_er = fid_er1632;
                    sm_win = sm_win1632;
                end
                df_3 = sprintf("%s/%s",chek,aw(j).name);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fid_i = fopen(df_3,'r');
                temp = fread(fid_i,'*double');
                fclose(fid_i);
                dta = temp(1:end-1);
                t_act = linspace(-1,(length(dta)/100)-1,length(dta));
                if (strfind(aw(j).name,'IDI'))
                    fsamp = 80;
                    t_len = length(dta)/fsamp;
                    l_nw = 100 * t_len;
                    x = linspace(0,t_len,length(dta));
                    x_nw = linspace(0,t_len,l_nw);
                    dta_nw = interp1(x,dta,x_nw);
                    dta = dta_nw';
                end
                if (isempty(strfind(aw(j).name,'Henvwin')))
                    if (strfind(aw(j).name,'Benvwin'))
                        fsamp = 20;
                        t_len = length(dta)/fsamp;
                        l_nw = 100 * t_len;
                        x = linspace(0,t_len,length(dta));
                        x_nw = linspace(0,t_len,l_nw);
                        dta_nw = interp1(x,dta,x_nw);
                        dta = dta_nw';
                    else
                        fsamp = 50;
                        t_len = length(dta)/fsamp;
                        l_nw = 100 * t_len;
                        x = linspace(0,t_len,length(dta));
                        x_nw = linspace(0,t_len,l_nw);
                        dta_nw = interp1(x,dta,x_nw);
                        dta = dta_nw';
                    end
                end
                bart = bartlett(sm_win);
                A_conv = conv(dta,bart,'same');%./100;
                A_conv = A_conv * max(dta)/max(A_conv);
                %pkg load signal
                [q,loc] = findpeaks(A_conv(1001:end));
                A_inv = 1./A_conv;
                [qd,locd] = findpeaks(A_inv(1001:end));
                loc = loc + 1000;
                locd = locd + 1000;
                x = 1:1:length(A_conv);
                jmp = zeros(length(loc),1);
                idx = zeros(length(loc),1);
                for k=1:length(loc)
                    sw = find(locd<loc(k));
                    jmp(k) = 0;
                    idx(k) = 0;
                    if (~isempty(sw))
                        idx(k) = locd(sw(end));
                        jmp(k) = A_conv(loc(k))/A_conv(idx(k));
                    end
                end
                rw = find(jmp > 4);
                cut_p = length(A_conv) - 1;
                if (~isempty(rw))
                    cut_p = idx(rw(1));
                end
                E_arr = single([A_conv(1:cut_p)',temp(end)]);
                temp2 = A_conv(1:cut_p)';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                t_idx = find(isnan(E_arr));
                E_arr(t_idx) = [];
                
                sz(j) = single(length(E_arr));
                
                q1 = split(aw(j).name,'_');
                stnm = q1(3);
                
                idx_st(j) = find(strcmp(st_db{1},stnm));
                st_lat(j) = st_db{2}(idx_st(j));
                st_lon(j) = st_db{3}(idx_st(j));
                
                [st_x,st_y,st_z] = geodetic2enu(st_lat(j),st_lon(j),0,...
                                                34,18,0,wgs84);
                r = single(sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + (ev_z-st_z)^2));
                
                t_end = (sz(j)-1-100)/100;
                %j
                %{
                    X0 = 1:1:(sz(j)-1);
                    Y1 = ones(sz(j)-1,1)*log(E_arr(105));
                    [X0,Y0] = intersections(X0,log(E_arr(1:end-1)),X0,Y1,'ROBUST');
                    w_len = round(max(X0));
                    if w_len > d_S_len
                        w_len = d_S_len;
                    end
                %}
                w_len = d_S_len;
                
                if (sz(j) > envl_min)
                    
                    E = ones((sz(j)-w_len),1,'single');
                    
                    E(1:(sz(j)-w_len)) = E_arr(w_len:end-1);
                    E = log(E);
                    
                    tm = single(linspace(0.01,t_end,...
                         sz(j)-w_len));
                    %r
                    
                    vel_fl = find(strcmp(velst{1},stnm));
                    if (vel_fl)
                        v_avg = single(velst{2}(vel_fl));
                    else
                        v_avg = single(3300);
                    end
                    
                    v_avg = v_avg/1000;
                     
                    %%{
                    
                    err_sum = g_func_wrpr_mtrx_full(g_rng,E,...
                        tm,v_avg,r,sz(j)-w_len);
                    
                    [~,g_idx] = mink(err_sum,1);
                    
                    if (g_idx < 100) && (g_idx > 1)
                        
                        g_rng2 = logspace(log10(g_rng(g_idx-1)),...
                            log10(g_rng(g_idx+1)),100);
                    
                        err_sum = g_func_wrpr_mtrx_full(g_rng2,E,...
                            tm,v_avg,r,sz(j)-w_len);
                        
                        [~,g_idx2] = mink(err_sum,1);
                        
                        g_st_w = g_rng2(g_idx2);
                    
                    else
                        
                        g_st_w = g_rng(g_idx);
                    end
                     %}
                    temp3 = g_func_env_tst(g_st_w,v_avg,r,sz(j)-w_len,...
                        0);
                    
                    G_f = single(temp3'.*10^(-9));
                    
                    G_f = log(G_f);
                    
                    X = [ones((sz(j)-w_len),1,'single') -tm'];
                    
                    b = X\(E-G_f);
                    
                    g_st_tmp(j) = g_st_w;
                    
                    E_calc = G_f - b(2).*tm' + b(1);
                    
                    [Ao,Bo] = intersections(tm',E_calc,tm',E);
                    
                    w_len = single(round(Ao(2)*100) + 100);
                    
                    E = ones((sz(j)-w_len),1,'single');
                    G_f = ones((sz(j)-w_len),1,'single');
                    
                    E(1) = mean(E_arr(101:w_len));
                    E(2:(sz(j)-w_len)) = E_arr((w_len+1):end-1);
                    E = log(E);
                    
                    S_idx = single(find(E_arr(1:end-1)==max(E_arr(1:end-1))));
                    
                    tm = [single((S_idx(1)-101)/100) ...
                        single(linspace(((w_len-101)/100+0.01),t_end,...
                         sz(j)-w_len-1))];
                    
                    temp3 = g_func_env_tst(g_st_w,v_avg,r,sz(j)-1-100,...
                        0);
                    
                    temp1 = temp3.*10^(-9);
                    
                    G_f(1) = mean(temp1(1:w_len-100));
                    G_f(2:(sz(j)-w_len)) = temp1((w_len-99):end);
                    
                    G_f = log(G_f);
                    
                    X = [ones((sz(j)-w_len),1,'single') -tm'];
                    
                    b = X\(E-G_f);
                    
                    g_st_tmp(j) = g_st_w;
                    
                    E_calc = G_f - b(2).*tm' + b(1);
                    
                    sum_tmp(j) = sum(abs((E-E_calc)))...
                                    /length(E) * 100;
                    err_tmp(j) = abs(exp(E(1)) - exp(E_calc(1)));
                    error_prom(j) = abs(err_tmp(j)/exp(E(1)))*100;
                    
                    if(b(2) > 0) %&& (error_prom(j) < 1)
                        
                        if (freq==1.5)
                            lgW_tmp12(j,:) = b;
                        elseif (freq==3)
                            lgW_tmp24(j,:) = b;
                        elseif (freq==6)
                            lgW_tmp48(j,:) = b;
                        elseif (freq==12)
                            lgW_tmp816(j,:) = b;
                        else
                            lgW_tmp1632(j,:) = b;
                        end
                        %%{
                        fprintf(fid_rp,"%f %f %f %f %f %f %f %f %f\n",...
                            ev_lat,ev_lon,ev_dep,st_db{2}(idx_st(j)),...
                            st_db{3}(idx_st(j)),v_avg,t_end,g_st_w,b(2));
                        %}
                        %%{
                        fprintf(fid_er,"%s %s %s %s\n",dr_lst(i).name,...
                            string(b(2)),string(error_prom(j)),...
                            string(sum_tmp(j)));
                        %}
                        %%{
                        frq_tmp = split(fw{5},'.');
                        frq_rng = frq_tmp{1} + " " + "Hz";
                        i,j
                        msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                        msft_sgn = sum(E-E_calc);
                        if (msft_sgn > 0)
                            E_calc = E_calc + msft;
                        else
                            E_calc = E_calc - msft;
                        end
                        t_fnl = linspace(-1,tm(end),sz(j)-1);
                        
                        h1 = figure(i+j*2); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                        scatter(tm(1),E(1),'r','LineWidth',2); plot(tm,E,'r','LineWidth',2); %legend(string(r));
                        hold on; scatter(tm(1),E_calc(1),'b','LineWidth',2); plot(tm,E_calc,'b','LineWidth',2); ...
                            text(20,-23,frq_rng,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
                            plot(t_fnl,log(E_arr(1:end-1)),'g'); xlabel('Time (in s)'); xlim([0 30]);...
                            ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                        
                        figure(i+j*20); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                            xlabel('Time (in s)');...
                            ylabel('Amplitude');...
                            text(20,max(dta)/2,frq_rng,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
                        plot(t_fnl,E_arr(1:end-1),'LineWidth',2); hold on;
                        plot(t_fnl,temp(1:(sz(j)-1)),'LineWidth',2); xlim([0 30]);...
                            legend('Smoothed S-wave env.','Actual S-wave env.'); hold off;
                        
                        %}
                    end
                end
            end
            %%{
            WR_tmp = lgW_tmp12(:,1);
            lgW_inv = mean(WR_tmp(WR_tmp>0));
            R_inv = lgW_tmp12(:,1) - lgW_inv;
            waitbar(i/length(dr_lst),f,sprintf('%12.0f:%s',i,dr_lst(i).name))
            
            pls_idx = find(WR_tmp > 0);
            
            if (length(pls_idx) > 2)
                o_tme = sprintf("%s/%s/%sH%sM%sS%s",evnt_info{3},...
                    evnt_info{4},evnt_info{5},evnt_info{6},evnt_info{7},...
                    evnt_info{8}(1:2));
                str = sprintf('%s ',string(idx_st(pls_idx)));
                param_set2 = [string(i) o_tme string(ev_lat) ...
                    string(ev_lon) string(ev_dep) string(lgW_inv) str];
                fprintf(fidW12,"%s %s %s %s %s %s %s\n",param_set2');
                param_set3 = [idx_st(pls_idx)' R_inv(pls_idx) ...
                                ones(length(pls_idx),1)*i];
                fprintf(fidR12,"%d %f %d\n",param_set3');
            end
            
            WR_tmp = lgW_tmp24(:,1);
            lgW_inv = mean(WR_tmp(WR_tmp>0));
            R_inv = lgW_tmp24(:,1) - lgW_inv;
            waitbar(i/length(dr_lst),f,sprintf('%12.0f:%s',i,dr_lst(i).name))
            
            pls_idx = find(WR_tmp > 0);
            
            if (length(pls_idx) > 2)
                o_tme = sprintf("%s/%s/%sH%sM%sS%s",evnt_info{3},...
                    evnt_info{4},evnt_info{5},evnt_info{6},evnt_info{7},...
                    evnt_info{8}(1:2));
                str = sprintf('%s ',string(idx_st(pls_idx)));
                param_set2 = [string(i) o_tme string(ev_lat) ...
                    string(ev_lon) string(ev_dep) string(lgW_inv) str];
                fprintf(fidW24,"%s %s %s %s %s %s %s\n",param_set2');
                param_set3 = [idx_st(pls_idx)' R_inv(pls_idx) ...
                                ones(length(pls_idx),1)*i];
                fprintf(fidR24,"%d %f %d\n",param_set3');
            end
            
            WR_tmp = lgW_tmp48(:,1);
            lgW_inv = mean(WR_tmp(WR_tmp>0));
            R_inv = lgW_tmp48(:,1) - lgW_inv;
            waitbar(i/length(dr_lst),f,sprintf('%12.0f:%s',i,dr_lst(i).name))
            
            pls_idx = find(WR_tmp > 0);
            
            if (length(pls_idx) > 2)
                o_tme = sprintf("%s/%s/%sH%sM%sS%s",evnt_info{3},...
                    evnt_info{4},evnt_info{5},evnt_info{6},evnt_info{7},...
                    evnt_info{8}(1:2));
                str = sprintf('%s ',string(idx_st(pls_idx)));
                param_set2 = [string(i) o_tme string(ev_lat) ...
                    string(ev_lon) string(ev_dep) string(lgW_inv) str];
                fprintf(fidW48,"%s %s %s %s %s %s %s\n",param_set2');
                param_set3 = [idx_st(pls_idx)' R_inv(pls_idx) ...
                                ones(length(pls_idx),1)*i];
                fprintf(fidR48,"%d %f %d\n",param_set3');
            end
            
            WR_tmp = lgW_tmp816(:,1);
            lgW_inv = mean(WR_tmp(WR_tmp>0));
            R_inv = lgW_tmp816(:,1) - lgW_inv;
            waitbar(i/length(dr_lst),f,sprintf('%12.0f:%s',i,dr_lst(i).name))
            
            pls_idx = find(WR_tmp > 0);
            
            if (length(pls_idx) > 2)
                o_tme = sprintf("%s/%s/%sH%sM%sS%s",evnt_info{3},...
                    evnt_info{4},evnt_info{5},evnt_info{6},evnt_info{7},...
                    evnt_info{8}(1:2));
                str = sprintf('%s ',string(idx_st(pls_idx)));
                param_set2 = [string(i) o_tme string(ev_lat) ...
                    string(ev_lon) string(ev_dep) string(lgW_inv) str];
                fprintf(fidW816,"%s %s %s %s %s %s %s\n",param_set2');
                param_set3 = [idx_st(pls_idx)' R_inv(pls_idx) ...
                                ones(length(pls_idx),1)*i];
                fprintf(fidR816,"%d %f %d\n",param_set3');
            end
            
            WR_tmp = lgW_tmp1632(:,1);
            lgW_inv = mean(WR_tmp(WR_tmp>0));
            R_inv = lgW_tmp1632(:,1) - lgW_inv;
            waitbar(i/length(dr_lst),f,sprintf('%12.0f:%s',i,dr_lst(i).name))
            
            pls_idx = find(WR_tmp > 0);
            
            if (length(pls_idx) > 2)
                o_tme = sprintf("%s/%s/%sH%sM%sS%s",evnt_info{3},...
                    evnt_info{4},evnt_info{5},evnt_info{6},evnt_info{7},...
                    evnt_info{8}(1:2));
                str = sprintf('%s ',string(idx_st(pls_idx)));
                param_set2 = [string(i) o_tme string(ev_lat) ...
                    string(ev_lon) string(ev_dep) string(lgW_inv) str];
                fprintf(fidW1632,"%s %s %s %s %s %s %s\n",param_set2');
                param_set3 = [idx_st(pls_idx)' R_inv(pls_idx) ...
                                ones(length(pls_idx),1)*i];
                fprintf(fidR1632,"%d %f %d\n",param_set3');
            end
            %}
            
            %%{
             g_st_tmp(60:end)
%             lgW_tmp12
%             lgW_tmp24
%             lgW_tmp48
%             lgW_tmp816
            sum_tmp(60:end)
            err_tmp(60:end)
            error_prom(60:end)
            %}
            
            %{
            geoscatter(ev_lat,ev_lon,*); hold on; geobasemap
            %}

        end
    end
    if getappdata(f,'canceling')
        break
    end
end
fclose(fidW12);
fclose(fidR12);

fclose(fid_rp12);
fclose(fid_er12);

fclose(fidW24);
fclose(fidR24);

fclose(fid_rp24);
fclose(fid_er24);

fclose(fidW48);
fclose(fidR48);

fclose(fid_rp48);
fclose(fid_er48);

fclose(fidW816);
fclose(fidR816);

fclose(fid_rp816);
fclose(fid_er816);

fclose(fidW1632);
fclose(fidR1632);

fclose(fid_rp1632);
fclose(fid_er1632);

toc

delete(f);