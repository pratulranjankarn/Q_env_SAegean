%% Assigning search directory for events and reading station coordinates
%% and average station velocities 
% clc;
% close all;
% clearvars;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

% Directory list for events
dr_lst = [%dir("/home/user/env_EGEL/husn_events/2010*"); ...
    %dir("/home/user/env_EGEL/husn_events/2011*"); ...
    %dir("/home/user/env_EGEL/husn_events/2010-09*"); ...
    %dir("/home/user/env_EGEL/husn_events/2010-10*"); ...
    %dir("/home/user/env_EGEL/husn_events/2010-11*"); ...
    %dir("/home/user/env_EGEL/husn_events/2010-12*");
    dir("/home/user/env_EGEL/husn_events/2010*"); ...
    dir("/home/user/env_EGEL/husn_events/2011*"); ...
    dir("/home/user/env_EGEL/husn_events/2012*"); ...
    dir("/home/user/env_EGEL/husn_events/2013*"); ...
    dir("/home/user/env_EGEL/husn_events/2014*"); ...
    dir("/home/user/env_EGEL/husn_events/2015*"); ...
    dir("/home/user/env_EGEL/husn_events/2016*"); ...
    dir("/home/user/env_EGEL/husn_events/2017*"); ... 
    dir("/home/user/env_EGEL/husn_events/2018*"); ...
    dir("/home/user/env_EGEL/husn_events/2019*"); ...
    dir("/home/user/env_EGEL/EGEL_events/200*"); ...
    dir("/home/user/env_EGEL/Cycnet_events/0*"); ...
    dir("/home/user/env_EGEL/kef_events/2*"); ...
    dir("/home/user/env_EGEL/cornoth_events/2*");];

% file for station coordinates fmt= 'Name' 'Latitude' 'Longitude'
fid = fopen('stlist_f.txt','r');
st_db = textscan(fid,'%s %f %f');
fclose(fid);

% file for station velocities fmt= 'Name' 'Avg. Velocity'
fid = fopen('st_vel.txt','r');
velst = textscan(fid,'%s %f');
fclose(fid);

%wkdir = '/home/user/env_EGEL/env_gauss_7.2,3.6,1.8,.90_pkdbl_wlnchk';

%mkdir(wkdir)

wgs84 = wgs84Ellipsoid('kilometers');

%% considering 100 Hz sampling frequency
% Smoothing window size
% sm_win12 = 800; sm_win24 = 400; sm_win48 = 200; sm_win816 = 100;sm_win1632 = 45;

% sm_win12 = 720; sm_win24 = 360; sm_win48 = 180; sm_win816 = 90;sm_win1632 = 45;

% sm_win12 = 600; sm_win24 = 300; sm_win48 = 150; sm_win816 = 75;sm_win1632 = 45;

sm_win12 = 500; sm_win24 = 250; sm_win48 = 125; sm_win816 = 62.5;sm_win1632 = 45;

% Minimum total length of envelope in samples considering 100 Hz
% sampling frequency
%envl_min = single(1700);

% Direct S-wave window length in samples
%d_S_len = single(900);


%% Computing for all frequency bands

% abs_pdir = sprintf('%s/ray_prmtrs12.txt',wkdir);
% fid_rp12 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_W_1-2.txt',wkdir);
% fidW12 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_R_1-2.txt',wkdir);
% fidR12 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_err_1-2.txt',wkdir);
% fid_er12 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/ray_prmtrs24.txt',wkdir);
% fid_rp24 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_W_2-4.txt',wkdir);
% fidW24 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_R_2-4.txt',wkdir);
% fidR24 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_err_2-4.txt',wkdir);
% fid_er24 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/ray_prmtrs48.txt',wkdir);
% fid_rp48 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_W_4-8.txt',wkdir);
% fidW48 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_R_4-8.txt',wkdir);
% fidR48 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_err_4-8.txt',wkdir);
% fid_er48 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/ray_prmtrs816.txt',wkdir);
% fid_rp816 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_W_8-16.txt',wkdir);
% fidW816 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_R_8-16.txt',wkdir);
% fidR816 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_err_8-16.txt',wkdir);
% fid_er816 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/ray_prmtrs1632.txt',wkdir);
% fid_rp1632 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_W_16-32.txt',wkdir);
% fidW1632 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_R_16-32.txt',wkdir);
% fidR1632 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_err_16-32.txt',wkdir);
% fid_er1632 = fopen(abs_pdir,'w');
% 
% abs_pdir = sprintf('%s/env_wn_1-2.txt',wkdir);
% fid_wn12 = fopen(abs_pdir,'w');
% abs_pdir = sprintf('%s/env_wn_2-4.txt',wkdir);
% fid_wn24 = fopen(abs_pdir,'w');
% abs_pdir = sprintf('%s/env_wn_4-8.txt',wkdir);
% fid_wn48 = fopen(abs_pdir,'w');
% abs_pdir = sprintf('%s/env_wn_8-16.txt',wkdir);
% fid_wn816 = fopen(abs_pdir,'w');
% abs_pdir = sprintf('%s/env_wn_16-32.txt',wkdir);
% fid_wn1632 = fopen(abs_pdir,'w');


% f = waitbar(0,'1','Name','Iterating through files for all freq',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% setappdata(f,'canceling',0);

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
            b_tmp = zeros(n_fls,1,'single');
            sum_tmp = zeros(n_fls,1,'single');
            err_tmp = zeros(n_fls,1,'single');
            error_prom = zeros(n_fls,1,'single');
            
            sz = zeros(n_fls,1,'single');
            
            st_lat = zeros(n_fls,1,'single');
            st_lon = zeros(n_fls,1,'single');
            idx_st = zeros(1,n_fls,'single');
            
            b = zeros(2,n_fls,'single');
            g_sti = 1e-4;
%             FigH = figure('Position', get(0, 'Screensize'));
            for j = 60:n_fls
                fw = split(aw(j).name,'_');
                if (strcmp(fw{5},'1-2.bin'))
                    freq = 1.5;
%                     fid_rp = fid_rp12;
%                     fid_er = fid_er12;
%                     fid_wn = fid_wn12;
                    sm_win = sm_win12;
                    plotid = 1;
                elseif (strcmp(fw{5},'2-4.bin'))
                    freq = 3;
%                     fid_rp = fid_rp24;
%                     fid_er = fid_er24;
%                     fid_wn = fid_wn24;
                    sm_win = sm_win24;
                    plotid = 2;
                elseif (strcmp(fw{5},'4-8.bin'))
                    freq = 6;
%                     fid_rp = fid_rp48;
%                     fid_er = fid_er48;
%                     fid_wn = fid_wn48;
                    sm_win = sm_win48;
                    plotid = 3;
                elseif (strcmp(fw{5},'8-16.bin'))
                    freq = 12;
%                     fid_rp = fid_rp816;
%                     fid_er = fid_er816;
%                     fid_wn = fid_wn816;
                    sm_win = sm_win816;
                    plotid = 4;
                else
                    freq = 24;
%                     fid_rp = fid_rp1632;
%                     fid_er = fid_er1632;
%                     fid_wn = fid_wn1632;
                    sm_win = sm_win1632;
                    plotid = 5;
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
                bart = bartlett(round(sm_win));
                A_conv = conv(dta,bart/sum(bart),'same');
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
                r = single(sqrt((ev_x-st_x)^2 + (ev_y-st_y)^2 + ...
                    (ev_z-st_z)^2));
                
                vel_fl = find(strcmp(velst{1},stnm));
                if (vel_fl)
                    v_avg = single(velst{2}(vel_fl));
                else
                    v_avg = single(3300);
                end
                
                v_avg = v_avg/1000;
                
                t_end = (sz(j)-1-100)/100;
                t_cons = single(r/v_avg);
                %j
                %{
                    X0 = 1:1:(sz(j)-1);
                    Y1 = ones(sz(j)-1,1)*log(E_arr(105));
                    [X0,Y0] = intersections(X0,log(E_arr(1:end-1)),X0,Y1,...
                                    'ROBUST');
                    w_len = round(max(X0));
                    if w_len > d_S_len
                        w_len = d_S_len;
                    end
                %}
                
                S_idx = single(find(E_arr(1:end-1)==max(E_arr(1:end-1))));
                
                w_len = min(S_idx)*2;
                
                envl_min = w_len + 800;
                
                %w_len = d_S_len;
                
                if (sz(j) > envl_min) && (min(S_idx) > 100)
                    
                    E = ones((sz(j)-w_len),1,'single');
                    G_f = ones((sz(j)-w_len),1,'single');
                    
                    E(1) = mean(E_arr(101:w_len));
                    E(2:(sz(j)-w_len)) = E_arr((w_len+1):end-1);
                    E = log(E);
                    
                    tm = [single((S_idx(1)-101)/100+t_cons) ...
                        single(linspace(((w_len-100)/100+0.01)+t_cons,...
                        t_end+t_cons,sz(j)-w_len-1))];
                    %r
                    
                    
                     
                    %%{
                    
                    err_sum = g_func_wrpr_mtrx(g_rng,w_len-100,E,...
                        tm,v_avg,r,sz(j)-1-99);
                    err_sum_tmp = err_sum;
                    
                    [~,g_idx] = mink(err_sum,1);
                    
                    if (g_idx < 100) && (g_idx > 1)
                        
                        g_rng2 = logspace(log10(g_rng(g_idx-1)),...
                            log10(g_rng(g_idx+1)),100);
                    
                        err_sum = g_func_wrpr_mtrx(g_rng2,w_len-100,E,...
                            tm,v_avg,r,sz(j)-1-99);
                        
                        [~,g_idx2] = mink(err_sum,1);
                        
                        g_st_w = g_rng2(g_idx2);
                    
                    else
                        
                        g_st_w = g_rng(g_idx);
                    end
                     %}
                    temp3 = g_func_env_tst(g_st_w,v_avg,r,sz(j)-1-100,...
                        0);
                    
                    temp1 = temp3.*10^(-9);
                    
                    G_f(1) = mean(temp1(1:w_len-100));
                    G_f(2:(sz(j)-w_len)) = temp1((w_len-99):end);
                    
                    G_f = log(G_f);
                    
                    X = [ones((sz(j)-w_len),1,'single') -tm'];
                    
                    b = X\(E-G_f);
                    
                    g_st_tmp(j) = g_st_w;
                    b_tmp(j) = b(2);
                    
                    E_calc = G_f - b(2).*tm' + b(1);
                    
                    sum_tmp(j) = sum(abs((E-E_calc)))...
                                    /length(E) * 100;
                    err_tmp(j) = abs(exp(E(1)) - exp(E_calc(1)));
                    error_prom(j) = abs(err_tmp(j)/exp(E(1)))*100;
                    
                    if(b(2) > 0) && (error_prom(j) < 1)
                        
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
                        %{
                        fprintf(fid_rp,"%f %f %f %f %f %f %f %f %f %f\n",...
                            ev_lat,ev_lon,ev_dep,st_db{2}(idx_st(j)),...
                            st_db{3}(idx_st(j)),v_avg,t_end,g_st_w,b(2),...
                            w_len);
                        %}
                        %{
                        fprintf(fid_er,"%s %s %s %s\n",dr_lst(i).name,...
                            string(b(2)),string(error_prom(j)),...
                            string(sum_tmp(j)));
                        %}
                        %{
                        fprintf(fid_wn,"%f\n",w_len);
                        
                        %}
                        %%{
                        frq_tmp = split(fw{5},'.');
                        frq_rng = frq_tmp{1} + " " + "Hz";
                        i,j,(i+j*2),(i+j*20), freq
                        
                        n_org = length(temp);
                        
                        t_fnl = linspace(t_cons+round(-1)+0.01,tm(end),...
                                            sz(j)-1);
                        t_end2 = (n_org-100-1)/100+t_cons;
                        t_fnl2 = linspace(t_cons+round(-1)+0.01,t_end2,...
                                            n_org-1);
                        
                        E_calc2 = zeros(1,n_org-100,'single');
                        
                        temp4 = g_func_env_tst(g_st_w,v_avg,r,n_org-1-100,0);
                        
                        temp3 = temp4.*10^(-9);
                        
                        tm2 = linspace(t_cons,t_end2,n_org-100);
                        E_calc2(2:end) = log(temp3) - b(2) .* tm2(2:end)...
                                            + b(1);
                        E_calc2(1) = -45;
                        
                        length(temp)
                        sz(j)
                       
                            %subplot(2,2,plotid); ...
                            if (rem(plotid,4) == 1)
%                                 h = axes('Position',[0.08 0.57 0.43 0.40]);
%                                 h = axes('Position',[0.08 0.74 0.21 0.18]);
%                                 h = axes('Position',[0.08 0.53 0.21 0.18]);
%                                 h = axes('Position',[0.08 0.32 0.21 0.18]);
                                h = axes('Position',[0.08 0.11 0.21 0.18]);
                                set(h,'Box','on',...
                                    'BoxStyle','full','LineWidth',3,...
                                    'XMinorTick','on','YMinorTick','on');
%                                 set(h,'XTickLabel',{});
                            elseif (rem(plotid,4) == 2)
%                                 h = axes('Position',[0.55 0.57 0.43 0.40]);
%                                 h = axes('Position',[0.31 0.74 0.21 0.18]);
%                                 h = axes('Position',[0.31 0.53 0.21 0.18]);
%                                 h = axes('Position',[0.31 0.32 0.21 0.18]);
                                h = axes('Position',[0.31 0.11 0.21 0.18]);
                                set(h,'YTickLabel',{});
%                                 set(h,'XTickLabel',{});
                                set(h,'Box','on',...
                                    'BoxStyle','full','LineWidth',3,...
                                    'XMinorTick','on','YMinorTick','on');
                            elseif (rem(plotid,4) == 3)
%                                 h = axes('Position',[0.08 0.12 0.43 0.40]);
%                                 h = axes('Position',[0.54 0.74 0.21 0.18]);
%                                 h = axes('Position',[0.54 0.53 0.21 0.18]);
%                                 h = axes('Position',[0.54 0.32 0.21 0.18]);
                                h = axes('Position',[0.54 0.11 0.21 0.18]);
                                set(h,'Box','on','BoxStyle','full',...
                                    'LineWidth',3,'XMinorTick','on'...
                                    ,'YMinorTick','on');
                                set(h,'YTickLabel',{});
%                                 set(h,'XTickLabel',{});
                            else
%                                 h = axes('Position',[0.55 0.12 0.43 0.4]);
%                                 h = axes('Position',[0.77 0.74 0.21 0.18]);
%                                 h = axes('Position',[0.77 0.53 0.21 0.18]);
%                                 h = axes('Position',[0.77 0.32 0.21 0.18]);
                                h = axes('Position',[0.77 0.11 0.21 0.18]);
                                set(h,'YTickLabel',{},'Box',...
                                    'on','BoxStyle','full','LineWidth',3,...
                                    'XMinorTick','on','YMinorTick','on');
                                set(h,'YTickLabel',{});
%                                 set(h,'XTickLabel',{});
                            end
                            set(h,'FontName','Helvetica','Fontsize',16,...
                                'FontWeight','bold'); hold on;...
%                             xlabel('Time (in s)');...
%                             ylabel('Log Amplitude');...
                            semilogy(t_fnl2,log(temp(1:end-1)),...
                                'LineWidth',2,'Color',[170 170 170]/255);...
                            hold on;...
                            blu_hsv = rgb2hsv([0 0 255]./255);
                            blu_hfstr_rgb = hsv2rgb([blu_hsv(1) blu_hsv(2)...
                                *0.5 blu_hsv(3)]);
                            semilogy(t_fnl2,log(A_conv),'LineWidth',4,...
                                'Color',blu_hfstr_rgb);%[111 168 220]/255*1);
                            semilogy(t_fnl2(w_len+1:cut_p),...
                                log(A_conv(w_len+1:cut_p)),'LineWidth',4,...
                                'Color',[0 0 255]/255*1);
                            rd_hsv = rgb2hsv([255 0 0]./255);
                            rd_hfstr_rgb = hsv2rgb([rd_hsv(1) rd_hsv(2)...
                                *0.4 rd_hsv(3)]);
                            %title(frq_rng);...
                            plot(tm2,E_calc2,'Color',rd_hfstr_rgb,...
                                'LineWidth',4); plot(tm(2:end),...
                                E_calc(2:end),'Color',[255 0 0]/255*1,...
                                'LineWidth',4); ...
                            plot(tm2(1:w_len-100),-40*ones(1,w_len-100),'Color',[176 194 74]/255*1,...
                                'LineWidth',7);
                            plot(tm(60:end),-40*ones(1,length(tm)-59),'Color',[176,194,74]/255*1,...
                                'LineWidth',7);
                            
%                             plot(tm2(1:w_len-100),E_calc(1)*ones(1,w_len-100),...
%                                 'Color',[0,0,255]/255*1,'LineWidth',3,...
%                             'LineStyle','-.');
                        
                            scatter(tm(1),E(1),'Marker','o',...
                                'MarkerEdgeColor',[128,0,128]/255,...
                                'LineWidth',2,'MarkerFaceColor',...
                                [128 0 128]/255,'SizeData',50); ...
                            xlim([0 84]);
                            ylim([-45 -10]);
                            Qsc_j = g_st_tmp(j) * v_avg/(2*pi*freq);
                            Qi_j = b_tmp(j)/(2*pi*freq);
                            Qscdsply = sprintf('Q_{sc}^{-1} =  %2.2e',Qsc_j);
                            Qidsply = sprintf('Q_{i}^{-1} =  %2.2e',Qi_j);
                            WRdsply = sprintf('log(WR) =  %2.2f',b(1));
%                             text(55,-15,Qscdsply,'FontName',...
%                             'Helvetica','Fontsize',18,'FontWeight','bold');
%                             text(55,-19,Qidsply,'FontName',...
%                             'Helvetica','Fontsize',18,'FontWeight','bold');
                            sbpltdsp = sprintf('W_{smooth} = %s s',...
                                string(sm_win/100));
                            text(27,-18,sbpltdsp,'FontName',...
                            'Helvetica','Fontsize',20,'FontWeight','bold');
                            
                            hold off;
                            
                            
                        %}
                        %{
                        figure(i+j*2); plot(err_sum_tmp);
                        
                        %}
%                         g_rng2(g_idx2)
%                         b(2)
%                         j
                    end
                end
            end
            %%{
            h4 = suplabel('Log(E) (Jm^{-3}Hz^{-1})','y');
            set(h4,'Position',[0.070 0.1300 0.7750 0.8150],'FontSize',20,...
                'FontWeight','bold');
            
            h5 = suplabel('Time since origin (in s)','x');
            set(h5,'Position',[0.135 0.1000 0.7750 0.8150],'FontSize',20,...
                'FontWeight','bold');
            %}
            %%{
            
            h5 = suplabel('1-2 Hz','x');
            set(h5,'Position',[0.1750 1.0000 0.0500 0.7150],'FontSize',18,...
                'FontWeight','bold');
            h5 = suplabel('2-4 Hz','x');
            set(h5,'Position',[0.1750 1.0000 0.5100 0.7150],'FontSize',18,...
                'FontWeight','bold');
            h5 = suplabel('4-8 Hz','x');
            set(h5,'Position',[0.1750 1.0000 0.9700 0.7150],'FontSize',18,...
                'FontWeight','bold');
            h5 = suplabel('8-16 Hz','x');
            set(h5,'Position',[0.1750 1.0000 1.4300 0.7150],'FontSize',18,...
                'FontWeight','bold');
            %}
            %waitbar(i/length(dr_lst),f,sprintf('%12.0f:%s',i,...
            %           dr_lst(i).name))
            %{
            WR_tmp = lgW_tmp12(:,1);
            lgW_inv = mean(WR_tmp(WR_tmp>0));
            R_inv = lgW_tmp12(:,1) - lgW_inv;
            
            
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
            
            %{
            g_st_tmp(60:end)
            b_tmp(60:end)
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
%     if getappdata(f,'canceling')
%         break
%     end
end
% fclose(fidW12);
% fclose(fidR12);
% 
% fclose(fid_rp12);
% fclose(fid_er12);
% 
% fclose(fidW24);
% fclose(fidR24);
% 
% fclose(fid_rp24);
% fclose(fid_er24);
% 
% fclose(fidW48);
% fclose(fidR48);
% 
% fclose(fid_rp48);
% fclose(fid_er48);
% 
% fclose(fidW816);
% fclose(fidR816);
% 
% fclose(fid_rp816);
% fclose(fid_er816);
% 
% fclose(fidW1632);
% fclose(fidR1632);
% 
% fclose(fid_rp1632);
% fclose(fid_er1632);
% 
% fclose(fid_wn12);
% 
% fclose(fid_wn24);
% 
% fclose(fid_wn48);
% 
% fclose(fid_wn816);
% 
% fclose(fid_wn1632);

toc

% delete(f);