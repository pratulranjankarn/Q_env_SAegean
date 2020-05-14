function temp1 = g_func_env(g_st,v_avg,r,sz,w_len,t_end,b)

% % clc;
% % close all;
% % clearvars;
% tic
% g_st = single(0.00053);
% v_avg = single(3.3373);
% r = single(10.830);
% sz = single(10000);
% w_len = single(1000);
% t_end = single(100);
% b = single(0.0053);

    temp1 = zeros(1,sz,'single');
    
    t_len = single(w_len/100);

    tcons = single(r/v_avg);
    
    
    tm = single(linspace(tcons,t_len+tcons,w_len+1));
    %length(tm)
    scndtrm_p5 = zeros(1,w_len+1,'single');
    
    fctr = exp(-v_avg*tm*g_st(1));
    frstrm = zeros(1,length(tm),'single');
    [~,qw] = min(tm-r/v_avg,[],2);
    frstrm(qw) = 1/(4*pi*r^2);
    scndtrm_p1 = exp(-1.5*(log(4*pi*v_avg)-log(3*g_st(1))));
    scndtrm_p2 = exp((-1.5)*log(tm));
    alpha = r./(v_avg.*tm);
    beta = (1 - alpha.^2);
    beta(qw) = 0;
    scndtrm_p3 = exp(log(beta)*(1/8));
    chek = v_avg*tm-r;
    scndtrm_p5(chek > 0) = 1;
    scndtrm_p5(chek == 0) = 0.5;
    scndtrm_p5(chek < 0) = 0;
    %scndtrm_p5 = heaviside(v_avg*tm-r);
    scndtrm_p4_var = exp(log(v_avg*tm*g_st(1)) + 0.75*log(beta));
    scndtrm_p4_var(qw) = 0;
    scndtrm_p4 = exp(scndtrm_p4_var + 0.5*log(1 + 2.026./scndtrm_p4_var));
    scndtrm_p4(qw) = 0;
    temp1(1:w_len+1) = fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5) .* exp(-b(1)*tm);
    
    %%{
    tm = single(linspace(t_len+0.01+tcons,t_end+tcons,(sz-w_len-1)));
    %length(tm)
    scndtrm_p5 = zeros(1,sz-w_len-1,'single');
    
    fctr = exp(-v_avg*tm*g_st(1));
    frstrm = single(0);
    scndtrm_p1 = (4*pi*v_avg/(3*g_st(1)))^(-1.5);
    scndtrm_p2 = exp((-1.5)*log(tm));
    scndtrm_p3 = exp(log(1 - r^2./(v_avg^2*tm.^2))*(1/8));
    scndtrm_p4_var = v_avg*tm*g_st(1).*(1-r^2./(v_avg^2*tm.^2)).^0.75;
    scndtrm_p4 = exp(scndtrm_p4_var).* sqrt(1 + 2.026./scndtrm_p4_var);
    chek = v_avg*tm-r;
    scndtrm_p5(chek > 0) = 1;
    scndtrm_p5(chek == 0) = 0.5;
    scndtrm_p5(chek < 0) = 0;
    %temp1(w_len+2:sz) = fctr .* (scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5) .* exp(-b(1)*tm);
    temp1(w_len+2:sz) = exp(log(fctr) + log(scndtrm_p1) + log(scndtrm_p2) + log(scndtrm_p3) + log(scndtrm_p4) + log(scndtrm_p5) -b(1)*tm);
    %}
    toc
    %plot(log(temp1(1:w_len+1)))
end