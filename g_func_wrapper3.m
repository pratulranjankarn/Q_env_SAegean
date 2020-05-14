function error_sum = g_func_wrapper3(g_st,w_len,E,tm,v_avg,r,sz)

    g_st = single(g_st);

    G_f = ones((sz-w_len),1,'single');
    
    t_len = single(sz/100);

    tcons = single(r/v_avg);
    
    
    tm_dm = tcons:0.01:(t_len+tcons);
    
    log_fctr = -v_avg*tm_dm*g_st(1);
    frstrm = 1/(4*pi*r^2);
    fctr_tmp1 = log_fctr(1) + log(frstrm);
    trm_1 = exp(fctr_tmp1);
    
    log_scndtrm_p1 = -1.5*(log(4*pi*v_avg)-log(3*g_st(1)));
    log_scndtrm_p2 = (-1.5)*log(tm_dm);
    log_scndtrm_p3 = log(1 - r^2./(v_avg^2*tm_dm.^2))*0.125;
    scndtrm_p4_var = exp(6*log_scndtrm_p3 + log(-log_fctr));
    log_scndtrm_p4 = scndtrm_p4_var + 0.5 * log(1 + 2.026./scndtrm_p4_var);
    fctr_tmp2 = log_scndtrm_p4 + log_scndtrm_p3 + ...
        log_scndtrm_p1 + log_scndtrm_p2 + log_fctr;
    trm_2 = exp(fctr_tmp2);
    trm_2(1) = trm_1;
    temp1 = trm_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_f1 = temp1.*10^(-9);

G_f(1) = sum(G_f1(1:w_len+1));
G_f(2:(sz-w_len)) = G_f1((w_len+2):end-1);

G_f = log(G_f);

X = [ones((sz-w_len),1,'single') -tm'];

% size(G_f1)
% size(G_f)
% size(X)
% size(E)


b = X\(E-G_f);

%{
%if (b(2) > 0.006) || (b(2) < 0.004)
    lb = [0 1e-5];
    ub = [inf 1e-1];
    options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
    options.Display = 'off';
    %options = optimoptions(options,'SubproblemAlgorithm','factorization');
    %options = optimoptions(options,'FunctionTolerance',1.000000e-10,'OptimalityTolerance',1e-34);
    [b,resnorm,residual,exitflag,output] = ...
        lsqlin(double(X),double((E-G_f)),[],[],[],[],lb,ub,[],options);
%end
%}
sum_tmp = ((E-G_f)-X*b).^2;
%E_calc = G_f+b(1,j)-b(2,j).*tm';
%t_fnl = linspace(-1.0,tm(end),sz-1);
%figure(j)
%scatter(tm(1),E(1)); hold on; plot(tm,E); %legend(string(r));
%hold on; scatter(tm(1),E_calc(1)); plot(tm,E_calc); plot(t_fnl,log(E_arr(1:end-1))); hold off;
error_sum = double((sum(sum_tmp(2:end))+w_len*sum_tmp(1)));
%%{
if (~isreal(b)) || (isempty(b)) %|| (~isreal(g_st))
    error_sum = inf;
end
%}
%clearvars('-except','error_sum');
end