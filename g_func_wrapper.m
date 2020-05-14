function error_sum = g_func_wrapper(aw,pw,g_st,w_len,E_arset,velst)
sum_tmp = cell(1,length(aw));
%sum_tmp=0;
g_st = single(g_st);
b = zeros(2,length(aw),'single');
    for j = 1:length(aw)
       % j
       E_arr = single(E_arset{j});
       sz = single(length(E_arr));
       if (sz > 1901)
       E = ones(1,(sz-w_len),'single');
       %G_f = ones(1,(sz-w_len),'single');
       tm = zeros((sz-w_len),1,'single');
       t_end = (sz-1-100)/100;
       %isreal(E_arr)
       E(1) = sum(E_arr(101:w_len))/(w_len-100)*1000;
       E(2:(sz-w_len)) = E_arr((w_len+1):end-1)*1000;
       r = E_arr(end)*1000;
       stnm = erase(pw{(j-1)*6+1},["outful_N","_Henvwin_4-8.txt","_"]);
       vel_fl = find(strcmp(velst{1},stnm));
       if (vel_fl)
           v_avg = single(velst{2}(vel_fl));
       else 
           v_avg = single(3300);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        temp1 = zeros(1,(sz-w_len),'single');
    
        t_len = single((w_len-100)/100);

        tcons = single(r/v_avg);
    
    
        tm_dm = single(linspace(tcons,t_len+tcons,(w_len-100)+1));
        %length(tm_dm)
    
        fctr = exp(-v_avg*tm_dm*g_st(j));
        frstrm = zeros(1,length(tm_dm));
        [~,qw] = min(tm_dm-r/v_avg,[],2);
        frstrm(qw) = 1/(4*pi*r^2);
        scndtrm_p1 = exp(-1.5*(log(4*pi*v_avg)-log(3*g_st(j))));
        scndtrm_p2 = exp((-1.5)*log(tm_dm));
        alpha = r./(v_avg.*tm_dm);
        beta = (1 - alpha.^2);
        beta(qw) = 0;
        scndtrm_p3 = exp(log(beta)*(1/8));
        chek = v_avg*tm_dm-r;
        if (chek > 0)
            scndtrm_p5 = 1;
        elseif (chek == 0)
            scndtrm_p5 = 1/2;
        else
            scndtrm_p5 = 0;
        end
        %scndtrm_p5 = heaviside(v_avg*tm_dm-r);
        scndtrm_p4_var = exp(log(v_avg*tm_dm*g_st(j)) + 0.75*log(beta));
        scndtrm_p4_var(qw) = 0;
        scndtrm_p4 = exp(scndtrm_p4_var + 0.5*log(1 + 2.026./scndtrm_p4_var));
        scndtrm_p4(qw) = 0;
        temp1(1) = sum(fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5))/((w_len-100)+1);
    
        tm_dm = linspace(t_len+0.01+tcons,t_end+tcons,((sz-w_len)-1));
        %length(tm_dm)
    
        fctr = exp(-v_avg*tm_dm*g_st(j));
        frstrm = 0;
        scndtrm_p1 = (4*pi*v_avg/(3*g_st(j)))^(-1.5);
        scndtrm_p2 = exp((-1.5)*log(tm_dm));
        scndtrm_p3 = exp(log(1 - r^2./(v_avg^2*tm_dm.^2))*(0.125));
        scndtrm_p4_var = exp(log(v_avg*tm_dm*g_st(j))+log(scndtrm_p3)*6);
        scndtrm_p4 = exp(scndtrm_p4_var).* sqrt(1 + 2.026./scndtrm_p4_var);
        chek = v_avg*tm_dm-r;
        if (chek > 0)
            scndtrm_p5 = 1;
        elseif (chek == 0)
            scndtrm_p5 = 1/2;
        else
            scndtrm_p5 = 0;
        end
        temp1(2:(sz-w_len)) = fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       G_f = log(temp1);
       E = log(E);
       S_idx = single(find(E_arr(1:w_len)==max(E_arr(1:w_len))));
       tm(1) = single((S_idx(1)-101)/100);
       %tm(1) = (median(find(abs(log(E_arr(1:end))-E(1)) < 0.01))-101)/100;
       tm(2:(sz-w_len)) = single(linspace(10.01,t_end,sz-1-w_len));
       %tm(E==0) = [];
       %E(E==0) = [];
       %G_f(G_f==0) = [];
       %length(E)
       %length(tm)
       %length(G_f)
       X = [ones(length(G_f),1,'single') -tm];
       b(:,j) = double(X)\(double(E-G_f))';
       sum_tmp{j} = double(sum(((E-G_f)'-X*b(:,j)).^2));
       %E_calc = G_f+b(1,j)-b(2,j).*tm';
       %t_fnl = linspace(-1.0,tm(end),sz-1);
       %figure(j)
       %scatter(tm(1),E(1)); hold on; plot(tm,E); %legend(string(r));
       %hold on; scatter(tm(1),E_calc(1)); plot(tm,E_calc); plot(t_fnl,log(E_arr(1:end-1))); hold off;
       end
    end
    qw = b(2,:)==0;
    b(:,qw) = [];
    %b
    error_sum = cell2mat(sum_tmp);
    if (~isreal(b)) || (isempty(b))
        error_sum = inf;
    end
    %clearvars('-except','sum_tmp');
end