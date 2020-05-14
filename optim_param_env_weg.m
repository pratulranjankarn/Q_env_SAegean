clc;
close all;
clearvars;


rw = [%dir("/home/user/env_EGEL/husn_events/2010*"); dir("/home/user/env_EGEL/husn_events/2011*");...
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

fid = fopen('stlist_f.txt','r');
st_db = textscan(fid,'%s %f %f');
fclose(fid);

fid = fopen('st_vel.txt','r');
velst = textscan(fid,'%s %f');
fclose(fid);

mw = struct2cell(rw);

clearvars('rw');

%%

fidWR = fopen('envweg_WR_1-2.txt','w');
fid_er = fopen('envweg_err_1-2.txt','w');

freq = 1.5;

%%%% for files with time length > 8s of coda considering 10s direct S window sze_lmt = 24700 bytes
sze_lmt = 24700;

w_len = single(1100);

for i=1:length(mw)
    chek = sprintf("%s/%s",mw{(i-1)*6+2},mw{(i-1)*6+1});
    chek
    df_1 = sprintf("%s/%s",chek,'outful*envwin_1-2.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        pw = struct2cell(aw);
        n_fl = length(aw);
        sze = zeros(n_fl,1);
        for k = 1:n_fl
            sze(k) = pw{(k-1)*6+4};
        end
        fl_idx = find(sze > sze_lmt);
        n_fls = length(fl_idx);
        if (n_fls > 0)
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            g_sti = 1e-4;
            
            df_2 = sprintf("%s/%s",chek,'env_param_weg_1-2.txt');
            fid = fopen(df_2,'w');
            Earset_cl = cell(1,n_fls);
            r = zeros(n_fls,1);
            v_avg = zeros(n_fls,1);
            S_idx = zeros(n_fls,1);
            sz = zeros(n_fls,1);
            tm_cl = cell(1,n_fls);
            stnm = cell(1,n_fls);
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                temp = load(df_3);
                t_idx = find(isnan(temp));
                temp(t_idx) = [];
                sze = length(temp);
                t_end = (sze-1-100)/100;
                sz(j) = sze;
                temp2 = [sum(temp(101:w_len))/(w_len-100)*1000 temp((w_len+1):end-1)'*1000];
                Earset_cl{j} = temp2;
                r(j) = temp(end)*1000;
                q1 = split(pw{(fl_idx(j)-1)*6+1},'_');
                stnm{j} = q1(3);
                vel_fl = find(strcmp(velst{1},stnm{j}));
                if (vel_fl)
                    v_avg(j) = single(velst{2}(vel_fl));
                else
                    v_avg(j) = single(3300);
                end
                idx_tmp = single(find(temp(1:w_len)==max(temp(1:w_len))));
                S_idx(j) = idx_tmp(1);
                tme = [single((S_idx(j)-101)/100) single(linspace(((w_len-100)/100+0.01),t_end,sze-1-w_len))];
                tm_cl{j} = tme;
            end
            
            Earset = cell2mat(Earset_cl);
            tm = cell2mat(tm_cl);
            
          
            fun = @(g_sti)g_func_wrapper_weg(g_sti,w_len,Earset,tm,v_avg,r,sz,n_fls);
            g_st = lsqnonlin(fun,g_sti,[],[],opts);
            %%{
            E = single(log(Earset));
            G_f_set = cell(1,n_fls);
            for j=1:n_fls
                v = v_avg(j);
                r_d = r(j);
                sz_d = sz(j);
                
                t_end = (sz_d-1-100)/100;
                
                temp1 = g_func_env_weg(g_st,v,r_d,(sz_d-w_len),(w_len-100),t_end);
                G_f_set{j} = log(temp1);
            end
            G_f = cell2mat(G_f_set);
            X = [ones(length(G_f),1,'single') -tm'];
            b = double(X)\(double(E-G_f))';
            if(~isnan(b))
                error_prom = 0;
                sum_tmp = sum(abs(((E-G_f)'-X*b)./E'));
                E_calc = G_f+b(1)-b(2).*tm;
                idx_i=1;
                idx_f=0;
                R = zeros(n_fls,1);
                for j=1:n_fls
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    error_prom = error_prom + abs((E_calc(idx_i) - E(idx_i))/E(idx_i));
                    idx_f = idx_f+sz(j)-w_len;
                    R(j) = exp(mean(E_calc(idx_i:idx_f)-E(idx_i:idx_f)));
                end     
                %if (error_prom < 1)
                param_set = [stnm' r/1000 string(ones(n_fls,1)*g_st)...
                    string(ones(n_fls,1)*b(2)) string(ones(n_fls,1)*exp(b(1)))...
                    string(R) ((sz-1-100)./100)];
                fprintf(fid,"%s %s %s %s %s %s %s\n",param_set');
                for j=1:n_fls
                    idx_st = find(strcmp(st_db{1},stnm{j}));
                    if(idx_st)
                        param_set2 = [i idx_st freq exp(b(1))];
                        fprintf(fidWR,"%f %f %f %f\n",param_set2');
                    end
                end
                fprintf(fid_er,"%s %s %s\n",mw{(i-1)*6+1},string(error_prom),string(sum_tmp));
                %end
                 %{
                 E_calc = G_f+b(1)-b(2).*tm;
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                idx_f=0;
                for j=1:n_fls
                    t_end = (sz(j)-1-100)/100;
                    t_fnl = linspace(-1.0,t_end,sz(j)-1);
                    
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    idx_f = idx_f+sz(j)-w_len;
                    
                    df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                    temp = load(df_3);
                    
                    h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                    scatter(tm(idx_i),E(idx_i),'r','LineWidth',2); plot(tm(idx_i:idx_f),E(idx_i:idx_f),'r','LineWidth',2); %legend(string(r));
                    hold on; scatter(tm(idx_i),E_calc(idx_i),'b','LineWidth',2); plot(tm(idx_i:idx_f),E_calc(idx_i:idx_f),'b','LineWidth',2); ...
                        plot(t_fnl,log(temp(1:end-1)*1000),'g');...
                        xlabel('Time (in s)');...
                        ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                    
                    figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                        xlabel('Time (in s)');...
                        ylabel('Amplitude');
                    plot(t_fnl,temp(1:end-1),'LineWidth',2); legend('Smoothed S-wave env.'); hold off;  
                    %hold off;
                end
                %}
            end
            fclose(fid);
        end
            %}
    end
end
fclose(fidWR);
fclose(fid_er);

%%

fidWR = fopen('envweg_WR_2-4.txt','w');
fid_er = fopen('envweg_err_2-4.txt','w');

freq = 3;

%%%% for files with time length > 8s of coda considering 10s direct S window sze_lmt = 24700 bytes
sze_lmt = 24700;

w_len = single(1100);

for i=1:length(mw)
    chek = sprintf("%s/%s",mw{(i-1)*6+2},mw{(i-1)*6+1});
    chek
    df_1 = sprintf("%s/%s",chek,'outful*envwin_2-4.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        pw = struct2cell(aw);
        n_fl = length(aw);
        sze = zeros(n_fl,1);
        for k = 1:n_fl
            sze(k) = pw{(k-1)*6+4};
        end
        fl_idx = find(sze > sze_lmt);
        n_fls = length(fl_idx);
        if (n_fls > 0)
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            g_sti = 1e-4;
            
            df_2 = sprintf("%s/%s",chek,'env_param_weg_2-4.txt');
            fid = fopen(df_2,'w');
            Earset_cl = cell(1,n_fls);
            r = zeros(n_fls,1);
            v_avg = zeros(n_fls,1);
            S_idx = zeros(n_fls,1);
            sz = zeros(n_fls,1);
            tm_cl = cell(1,n_fls);
            stnm = cell(1,n_fls);
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                temp = load(df_3);
                t_idx = find(isnan(temp));
                temp(t_idx) = [];
                sze = length(temp);
                t_end = (sze-1-100)/100;
                sz(j) = sze;
                temp2 = [sum(temp(101:w_len))/(w_len-100)*1000 temp((w_len+1):end-1)'*1000];
                Earset_cl{j} = temp2;
                r(j) = temp(end)*1000;
                q1 = split(pw{(fl_idx(j)-1)*6+1},'_');
                stnm{j} = q1(3);
                vel_fl = find(strcmp(velst{1},stnm{j}));
                if (vel_fl)
                    v_avg(j) = single(velst{2}(vel_fl));
                else
                    v_avg(j) = single(3300);
                end
                idx_tmp = single(find(temp(1:w_len)==max(temp(1:w_len))));
                S_idx(j) = idx_tmp(1);
                tme = [single((S_idx(j)-101)/100) single(linspace(((w_len-100)/100+0.01),t_end,sze-1-w_len))];
                tm_cl{j} = tme;
            end
            
            Earset = cell2mat(Earset_cl);
            tm = cell2mat(tm_cl);
            
          
            fun = @(g_sti)g_func_wrapper_weg(g_sti,w_len,Earset,tm,v_avg,r,sz,n_fls);
            g_st = lsqnonlin(fun,g_sti,[],[],opts);
            %%{
            E = single(log(Earset));
            G_f_set = cell(1,n_fls);
            for j=1:n_fls
                v = v_avg(j);
                r_d = r(j);
                sz_d = sz(j);
                
                t_end = (sz_d-1-100)/100;
                
                temp1 = g_func_env_weg(g_st,v,r_d,(sz_d-w_len),(w_len-100),t_end);
                G_f_set{j} = log(temp1);
            end
            G_f = cell2mat(G_f_set);
            X = [ones(length(G_f),1,'single') -tm'];
            b = double(X)\(double(E-G_f))';
            if(~isnan(b))
                error_prom = 0;
                sum_tmp = sum(abs(((E-G_f)'-X*b)./E'));
                E_calc = G_f+b(1)-b(2).*tm;
                idx_i=1;
                idx_f=0;
                R = zeros(n_fls,1);
                for j=1:n_fls
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    error_prom = error_prom + abs((E_calc(idx_i) - E(idx_i))/E(idx_i));
                    idx_f = idx_f+sz(j)-w_len;
                    R(j) = exp(mean(E_calc(idx_i:idx_f)-E(idx_i:idx_f)));
                end
                
                %if (error_prom < 1)
                param_set = [stnm' r/1000 string(ones(n_fls,1)*g_st)...
                    string(ones(n_fls,1)*b(2)) string(ones(n_fls,1)*exp(b(1)))...
                    string(R) ((sz-1-100)./100)];
                fprintf(fid,"%s %s %s %s %s %s %s\n",param_set');
                for j=1:n_fls
                    idx_st = find(strcmp(st_db{1},stnm{j}));
                    if(idx_st)
                        param_set2 = [i idx_st freq exp(b(1))];
                        fprintf(fidWR,"%f %f %f %f\n",param_set2');
                    end
                end
                fprintf(fid_er,"%s %s %s\n",mw{(i-1)*6+1},string(error_prom),string(sum_tmp));
                %end
                %{
                E_calc = G_f+b(1)-b(2).*tm;
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                idx_f=0;
                for j=1:n_fls
                    t_end = (sz(j)-1-100)/100;
                    t_fnl = linspace(-1.0,t_end,sz(j)-1);
                    
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    idx_f = idx_f+sz(j)-w_len;
                    
                    df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                    temp = load(df_3);
                    
                    h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                    scatter(tm(idx_i),E(idx_i),'r','LineWidth',2); plot(tm(idx_i:idx_f),E(idx_i:idx_f),'r','LineWidth',2); %legend(string(r));
                    hold on; scatter(tm(idx_i),E_calc(idx_i),'b','LineWidth',2); plot(tm(idx_i:idx_f),E_calc(idx_i:idx_f),'b','LineWidth',2); ...
                        plot(t_fnl,log(temp(1:end-1)*1000),'g');...
                        xlabel('Time (in s)');...
                        ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                    
                    figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                        xlabel('Time (in s)');...
                        ylabel('Amplitude');
                    plot(t_fnl,temp(1:end-1),'LineWidth',2); legend('Smoothed S-wave env.'); hold off;  
                    %hold off;
                end
                %}
                %hold off;
            end
            fclose(fid);
        end
            %}
    end
end
fclose(fidWR);
fclose(fid_er);

%%

fidWR = fopen('envweg_WR_4-8.txt','w');
fid_er = fopen('envweg_err_4-8.txt','w');

freq = 6;

%%%% for files with time length > 8s of coda considering 10s direct S window sze_lmt = 24700 bytes
sze_lmt = 24700;

w_len = single(1100);

for i=1:length(mw)
    chek = sprintf("%s/%s",mw{(i-1)*6+2},mw{(i-1)*6+1});
    chek
    df_1 = sprintf("%s/%s",chek,'outful*envwin_4-8.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        pw = struct2cell(aw);
        n_fl = length(aw);
        sze = zeros(n_fl,1);
        for k = 1:n_fl
            sze(k) = pw{(k-1)*6+4};
        end
        fl_idx = find(sze > sze_lmt);
        n_fls = length(fl_idx);
        if (n_fls > 0)
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            g_sti = 1e-4;
            
            df_2 = sprintf("%s/%s",chek,'env_param_weg_4-8.txt');
            fid = fopen(df_2,'w');
            Earset_cl = cell(1,n_fls);
            r = zeros(n_fls,1);
            v_avg = zeros(n_fls,1);
            S_idx = zeros(n_fls,1);
            sz = zeros(n_fls,1);
            tm_cl = cell(1,n_fls);
            stnm = cell(1,n_fls);
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                temp = load(df_3);
                t_idx = find(isnan(temp));
                temp(t_idx) = [];
                sze = length(temp);
                t_end = (sze-1-100)/100;
                sz(j) = sze;
                temp2 = [sum(temp(101:w_len))/(w_len-100)*1000 temp((w_len+1):end-1)'*1000];
                Earset_cl{j} = temp2;
                r(j) = temp(end)*1000;
                q1 = split(pw{(fl_idx(j)-1)*6+1},'_');
                stnm{j} = q1(3);
                vel_fl = find(strcmp(velst{1},stnm{j}));
                if (vel_fl)
                    v_avg(j) = single(velst{2}(vel_fl));
                else
                    v_avg(j) = single(3300);
                end
                idx_tmp = single(find(temp(1:w_len)==max(temp(1:w_len))));
                S_idx(j) = idx_tmp(1);
                tme = [single((S_idx(j)-101)/100) single(linspace(((w_len-100)/100+0.01),t_end,sze-1-w_len))];
                tm_cl{j} = tme;
            end
            
            Earset = cell2mat(Earset_cl);
            tm = cell2mat(tm_cl);
            
          
            fun = @(g_sti)g_func_wrapper_weg(g_sti,w_len,Earset,tm,v_avg,r,sz,n_fls);
            g_st = lsqnonlin(fun,g_sti,[],[],opts);
            %%{
            E = single(log(Earset));
            G_f_set = cell(1,n_fls);
            for j=1:n_fls
                v = v_avg(j);
                r_d = r(j);
                sz_d = sz(j);
                
                t_end = (sz_d-1-100)/100;
                
                temp1 = g_func_env_weg(g_st,v,r_d,(sz_d-w_len),(w_len-100),t_end);
                G_f_set{j} = log(temp1);
            end
            G_f = cell2mat(G_f_set);
            X = [ones(length(G_f),1,'single') -tm'];
            b = double(X)\(double(E-G_f))';
            if(~isnan(b))
                error_prom = 0;
                sum_tmp = sum(abs(((E-G_f)'-X*b)./E'));
                E_calc = G_f+b(1)-b(2).*tm;
                idx_i=1;
                idx_f=0;
                R = zeros(n_fls,1);
                for j=1:n_fls
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    error_prom = error_prom + abs((E_calc(idx_i) - E(idx_i))/E(idx_i));
                    idx_f = idx_f+sz(j)-w_len;
                    R(j) = exp(mean(E_calc(idx_i:idx_f)-E(idx_i:idx_f)));
                end
                
                %if (error_prom < 1)
                param_set = [stnm' r/1000 string(ones(n_fls,1)*g_st)...
                    string(ones(n_fls,1)*b(2)) string(ones(n_fls,1)*exp(b(1)))...
                    string(R) ((sz-1-100)./100)];
                fprintf(fid,"%s %s %s %s %s %s %s\n",param_set');
                for j=1:n_fls
                    idx_st = find(strcmp(st_db{1},stnm{j}));
                    if(idx_st)
                        param_set2 = [i idx_st freq exp(b(1))];
                        fprintf(fidWR,"%f %f %f %f\n",param_set2');
                    end
                end
                fprintf(fid_er,"%s %s %s\n",mw{(i-1)*6+1},string(error_prom),string(sum_tmp));
                %end
                 %{
                 E_calc = G_f+b(1)-b(2).*tm;
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                idx_f=0;
                for j=1:n_fls
                    t_end = (sz(j)-1-100)/100;
                    t_fnl = linspace(-1.0,t_end,sz(j)-1);
                    
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    idx_f = idx_f+sz(j)-w_len;
                    
                    df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                    temp = load(df_3);
                    
                    h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                    scatter(tm(idx_i),E(idx_i),'r','LineWidth',2); plot(tm(idx_i:idx_f),E(idx_i:idx_f),'r','LineWidth',2); %legend(string(r));
                    hold on; scatter(tm(idx_i),E_calc(idx_i),'b','LineWidth',2); plot(tm(idx_i:idx_f),E_calc(idx_i:idx_f),'b','LineWidth',2); ...
                        plot(t_fnl,log(temp(1:end-1)*1000),'g');...
                        xlabel('Time (in s)');...
                        ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                    
                    figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                        xlabel('Time (in s)');...
                        ylabel('Amplitude');
                    plot(t_fnl,temp(1:end-1),'LineWidth',2); legend('Smoothed S-wave env.'); hold off;  
                    %hold off;
                end
                %}
            end
            fclose(fid);
        end
        %}
    end
end
fclose(fidWR);
fclose(fid_er);

%%
fidWR = fopen('envweg_WR_8-16.txt','w');
fid_er = fopen('envweg_err_8-16.txt','w');

freq = 12;

%%%% for files with time length > 8s of coda considering 10s direct S window sze_lmt = 24700 bytes
sze_lmt = 24700;

w_len = single(1100);

for i=1:length(mw)
    chek = sprintf("%s/%s",mw{(i-1)*6+2},mw{(i-1)*6+1});
    chek
    df_1 = sprintf("%s/%s",chek,'outful*envwin_8-16.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        pw = struct2cell(aw);
        n_fl = length(aw);
        sze = zeros(n_fl,1);
        for k = 1:n_fl
            sze(k) = pw{(k-1)*6+4};
        end
        fl_idx = find(sze > sze_lmt);
        n_fls = length(fl_idx);
        if (n_fls > 0)
            
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            g_sti = 1e-4;
            
            df_2 = sprintf("%s/%s",chek,'env_param_weg_8-16.txt');
            fid = fopen(df_2,'w');
            Earset_cl = cell(1,n_fls);
            r = zeros(n_fls,1);
            v_avg = zeros(n_fls,1);
            S_idx = zeros(n_fls,1);
            sz = zeros(n_fls,1);
            tm_cl = cell(1,n_fls);
            stnm = cell(1,n_fls);
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                temp = load(df_3);
                t_idx = find(isnan(temp));
                temp(t_idx) = [];
                sze = length(temp);
                t_end = (sze-1-100)/100;
                sz(j) = sze;
                temp2 = [sum(temp(101:w_len))/(w_len-100)*1000 temp((w_len+1):end-1)'*1000];
                Earset_cl{j} = temp2;
                r(j) = temp(end)*1000;
                q1 = split(pw{(fl_idx(j)-1)*6+1},'_');
                stnm{j} = q1(3);
                vel_fl = find(strcmp(velst{1},stnm{j}));
                if (vel_fl)
                    v_avg(j) = single(velst{2}(vel_fl));
                else
                    v_avg(j) = single(3300);
                end
                idx_tmp = single(find(temp(1:w_len)==max(temp(1:w_len))));
                S_idx(j) = idx_tmp(1);
                tme = [single((S_idx(j)-101)/100) single(linspace(((w_len-100)/100+0.01),t_end,sze-1-w_len))];
                tm_cl{j} = tme;
            end
            
            Earset = cell2mat(Earset_cl);
            tm = cell2mat(tm_cl);
            
          
            fun = @(g_sti)g_func_wrapper_weg(g_sti,w_len,Earset,tm,v_avg,r,sz,n_fls);
            g_st = lsqnonlin(fun,g_sti,[],[],opts);
            %%{
            E = single(log(Earset));
            G_f_set = cell(1,n_fls);
            for j=1:n_fls
                v = v_avg(j);
                r_d = r(j);
                sz_d = sz(j);
                
                t_end = (sz_d-1-100)/100;
                
                temp1 = g_func_env_weg(g_st,v,r_d,(sz_d-w_len),(w_len-100),t_end);
                G_f_set{j} = log(temp1);
            end
            G_f = cell2mat(G_f_set);
            X = [ones(length(G_f),1,'single') -tm'];
            b = double(X)\(double(E-G_f))';
            if(~isnan(b))
                error_prom = 0;
                sum_tmp = sum(abs(((E-G_f)'-X*b)./E'));
                E_calc = G_f+b(1)-b(2).*tm;
                idx_i=1;
                idx_f=0;
                R = zeros(n_fls,1);
                for j=1:n_fls
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    error_prom = error_prom + abs((E_calc(idx_i) - E(idx_i))/E(idx_i));
                    idx_f = idx_f+sz(j)-w_len;
                    R(j) = exp(mean(E_calc(idx_i:idx_f)-E(idx_i:idx_f)));
                end
                error_prom = abs(E_calc(1) - E(1));
                %if (error_prom < 1)
                param_set = [stnm' r/1000 string(ones(n_fls,1)*g_st)...
                    string(ones(n_fls,1)*b(2)) string(ones(n_fls,1)*exp(b(1)))...
                    string(R) ((sz-1-100)./100)];
                fprintf(fid,"%s %s %s %s %s %s %s\n",param_set');
                for j=1:n_fls
                    idx_st = find(strcmp(st_db{1},stnm{j}));
                    if(idx_st)
                        param_set2 = [i idx_st freq exp(b(1))];
                        fprintf(fidWR,"%f %f %f %f\n",param_set2');
                    end
                end
                fprintf(fid_er,"%s %s %s\n",mw{(i-1)*6+1},string(error_prom),string(sum_tmp));
                %end
                %{
                 E_calc = G_f+b(1)-b(2).*tm;
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                idx_f=0;
                for j=1:n_fls
                    t_end = (sz(j)-1-100)/100;
                    t_fnl = linspace(-1.0,t_end,sz(j)-1);
                    
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    idx_f = idx_f+sz(j)-w_len;
                    
                    df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                    temp = load(df_3);
                    
                    h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                    scatter(tm(idx_i),E(idx_i),'r','LineWidth',2); plot(tm(idx_i:idx_f),E(idx_i:idx_f),'r','LineWidth',2); %legend(string(r));
                    hold on; scatter(tm(idx_i),E_calc(idx_i),'b','LineWidth',2); plot(tm(idx_i:idx_f),E_calc(idx_i:idx_f),'b','LineWidth',2); ...
                        plot(t_fnl,log(temp(1:end-1)*1000),'g');...
                        xlabel('Time (in s)');...
                        ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                    
                    figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                        xlabel('Time (in s)');...
                        ylabel('Amplitude');
                    plot(t_fnl,temp(1:end-1),'LineWidth',2); legend('Smoothed S-wave env.'); hold off;  
                    %hold off;
                end
                %}
            end
            fclose(fid);
        end
            %}
    end
end
fclose(fidWR);
fclose(fid_er);
%%


fidWR = fopen('envweg_WR_16-32.txt','w');
fid_er = fopen('envweg_err_16-32.txt','w');

freq = 24;

%%%% for files with time length > 8s of coda considering 10s direct S window sze_lmt = 24700 bytes
sze_lmt = 24700;

w_len = single(1100);

for i=1:length(mw)
    chek = sprintf("%s/%s",mw{(i-1)*6+2},mw{(i-1)*6+1});
    chek
    df_1 = sprintf("%s/%s",chek,'outful*envwin_16-32.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        pw = struct2cell(aw);
        n_fl = length(aw);
        sze = zeros(n_fl,1);
        for k = 1:n_fl
            sze(k) = pw{(k-1)*6+4};
        end
        fl_idx = find(sze > sze_lmt);
        n_fls = length(fl_idx);
        if (n_fls > 0)
            
            opts = optimoptions('lsqnonlin');
            opts.Algorithm = 'levenberg-marquardt';
            options.MaxIterations = 5000;
            g_sti = 1e-4;
            
            df_2 = sprintf("%s/%s",chek,'env_param_weg_16-32.txt');
            fid = fopen(df_2,'w');
            Earset_cl = cell(1,n_fls);
            r = zeros(n_fls,1);
            v_avg = zeros(n_fls,1);
            S_idx = zeros(n_fls,1);
            sz = zeros(n_fls,1);
            tm_cl = cell(1,n_fls);
            stnm = cell(1,n_fls);
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                temp = load(df_3);
                t_idx = find(isnan(temp));
                temp(t_idx) = [];
                sze = length(temp);
                t_end = (sze-1-100)/100;
                sz(j) = sze;
                temp2 = [sum(temp(101:w_len))/(w_len-100)*1000 temp((w_len+1):end-1)'*1000];
                Earset_cl{j} = temp2;
                r(j) = temp(end)*1000;
                q1 = split(pw{(fl_idx(j)-1)*6+1},'_');
                stnm{j} = q1(3);
                vel_fl = find(strcmp(velst{1},stnm{j}));
                if (vel_fl)
                    v_avg(j) = single(velst{2}(vel_fl));
                else
                    v_avg(j) = single(3300);
                end
                idx_tmp = single(find(temp(1:w_len)==max(temp(1:w_len))));
                S_idx(j) = idx_tmp(1);
                tme = [single((S_idx(j)-101)/100) single(linspace(((w_len-100)/100+0.01),t_end,sze-1-w_len))];
                tm_cl{j} = tme;
            end
            
            Earset = cell2mat(Earset_cl);
            tm = cell2mat(tm_cl);
            
          
            fun = @(g_sti)g_func_wrapper_weg(g_sti,w_len,Earset,tm,v_avg,r,sz,n_fls);
            g_st = lsqnonlin(fun,g_sti,[],[],opts);
            %%{
            E = single(log(Earset));
            G_f_set = cell(1,n_fls);
            for j=1:n_fls
                v = v_avg(j);
                r_d = r(j);
                sz_d = sz(j);
                
                t_end = (sz_d-1-100)/100;
                
                temp1 = g_func_env_weg(g_st,v,r_d,(sz_d-w_len),(w_len-100),t_end);
                G_f_set{j} = log(temp1);
            end
            G_f = cell2mat(G_f_set);
            X = [ones(length(G_f),1,'single') -tm'];
            b = double(X)\(double(E-G_f))';
            if(~isnan(b))
                error_prom = 0;
                sum_tmp = sum(abs(((E-G_f)'-X*b)./E'));
                E_calc = G_f+b(1)-b(2).*tm;
                idx_i=1;
                idx_f=0;
                R = zeros(n_fls,1);
                for j=1:n_fls
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    error_prom = error_prom + abs((E_calc(idx_i) - E(idx_i))/E(idx_i));
                    idx_f = idx_f+sz(j)-w_len;
                    R(j) = exp(mean(E_calc(idx_i:idx_f)-E(idx_i:idx_f)));
                end
                
                %if (error_prom < 1)
                param_set = [stnm' r/1000 string(ones(n_fls,1)*g_st)...
                    string(ones(n_fls,1)*b(2)) string(ones(n_fls,1)*exp(b(1)))...
                    string(R) ((sz-1-100)./100)];
                fprintf(fid,"%s %s %s %s %s %s %s\n",param_set');
                for j=1:n_fls
                    idx_st = find(strcmp(st_db{1},stnm{j}));
                    if(idx_st)
                        param_set2 = [i idx_st freq exp(b(1))];
                        fprintf(fidWR,"%f %f %f %f\n",param_set2');
                    end
                end
                fprintf(fid_er,"%s %s %s\n",mw{(i-1)*6+1},string(error_prom),string(sum_tmp));
                %end
                %{
                E_calc = G_f+b(1)-b(2).*tm;
                msft = sqrt(sum((E - E_calc).^2)/length(E_calc));
                msft_sgn = sum(E-E_calc);
                if (msft_sgn > 0)
                    E_calc = E_calc + msft;
                else
                    E_calc = E_calc - msft;
                end
                idx_f=0;
                for j=1:n_fls
                    t_end = (sz(j)-1-100)/100;
                    t_fnl = linspace(-1.0,t_end,sz(j)-1);
                    
                    if (j==1)
                        idx_i = 1;
                    else
                        idx_i = idx_i+sz(j-1)-w_len;
                    end
                    idx_f = idx_f+sz(j)-w_len;
                    
                    df_3 = sprintf("%s/%s",chek,pw{(fl_idx(j)-1)*6+1});
                    temp = load(df_3);
                    
                    h1 = figure(j*15); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;
                    scatter(tm(idx_i),E(idx_i),'r','LineWidth',2); plot(tm(idx_i:idx_f),E(idx_i:idx_f),'r','LineWidth',2); %legend(string(r));
                    hold on; scatter(tm(idx_i),E_calc(idx_i),'b','LineWidth',2); plot(tm(idx_i:idx_f),E_calc(idx_i:idx_f),'b','LineWidth',2); ...
                        plot(t_fnl,log(temp(1:end-1)*1000),'g');...
                        xlabel('Time (in s)');...
                        ylabel('log_{10}(E(r,t))'); legend('S-Winavg_{obs}','S-Coda_{obs}','S-Winavg_{calc}','S-Coda_{calc}'); hold off;
                    
                    figure(j*150); set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold'); hold on;...
                        xlabel('Time (in s)');...
                        ylabel('Amplitude');
                    plot(t_fnl,temp(1:end-1),'LineWidth',2); legend('Smoothed S-wave env.'); hold off;   
                    %hold off;
                end
                %}
            end
            fclose(fid);
        end
            %}
    end
end
fclose(fidWR);
fclose(fid_er);

