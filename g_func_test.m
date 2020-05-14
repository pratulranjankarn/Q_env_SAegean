%function temp1 = g_func_env_tst(g_st,v_avg,r,sz,b)
clc;
close all;
clearvars;
tic

tic
g_st = single(logspace(-7,1,100));
%g_st = g_st';
v_avg = single(3.3373);
r = single(40.830);
sz = single(10000);
w_len = single(1000);
t_end = single(100);
b = single(0.028)*0;

    %temp1 = zeros(1,sz,'single');
    
    t_len = single(sz/100);

    tcons = single(r/v_avg);
    
    
    tm = (tcons:0.01:(t_len+tcons-0.01))';
    %length(tm)
    %scndtrm_p5 = zeros(sz,1,'single');
    
    %trm_1 = zeros(sz,1,'single');
    
    log_fctr = -v_avg*tm*g_st;
    frstrm = 1/(4*pi*r^2);
    fctr_tmp1 = log_fctr(1) + log(frstrm) - b(1) * tm(1);
    trm_1 = exp(fctr_tmp1);
    
    log_scndtrm_p1 = -1.5*(log(4*pi*v_avg)-log(3*g_st));
    log_scndtrm_p2 = (-1.5)*log(tm);
    log_scndtrm_p3 = log(1 - r^2./(v_avg^2*tm.^2))*0.125;
    log_cmb_123 = log_scndtrm_p2 + log_scndtrm_p3 + log_scndtrm_p1;
    scndtrm_p4_var = (exp(6*log_scndtrm_p3)*v_avg.*tm)* g_st;
    log_scndtrm_p4 = scndtrm_p4_var + 0.5 * log(bsxfun(@plus,1,...
        bsxfun(@rdivide,2.026,scndtrm_p4_var)));
    fctr_tmp2 = log_scndtrm_p4 + log_cmb_123 + log_fctr;
    trm_2 = exp(fctr_tmp2- b(1)*tm);
    trm_2(1) = trm_1;
    temp1 = trm_2;
    
    figure(1);
    
    for i=1:100
    plot(log(temp1(10:10000,i)));
    hold on;
    end
    hold off;
    
    figure(3);
    
    plot(temp1(500,:));
    
    g_tmp = g_st;
    
    clearvars -except g_tmp
    
    g_st = single(logspace(log10(g_tmp(70)),log10(g_tmp(80)),100));
    %g_st = g_st';
    v_avg = single(3.3373);
    r = single(40.830);
    sz = single(10000);
    w_len = single(1000);
    t_end = single(100);
    b = single(0.028)*0;

    %temp1 = zeros(1,sz,'single');
    
    t_len = single(sz/100);

    tcons = single(r/v_avg);
    
    
    tm = (tcons:0.01:(t_len+tcons-0.01))';
    %length(tm)
    %scndtrm_p5 = zeros(sz,1,'single');
    
    %trm_1 = zeros(sz,1,'single');
    
    log_fctr = -v_avg*tm*g_st;
    frstrm = 1/(4*pi*r^2);
    fctr_tmp1 = log_fctr(1) + log(frstrm) - b(1) * tm(1);
    trm_1 = exp(fctr_tmp1);
    
    log_scndtrm_p1 = -1.5*(log(4*pi*v_avg)-log(3*g_st));
    log_scndtrm_p2 = (-1.5)*log(tm);
    log_scndtrm_p3 = log(1 - r^2./(v_avg^2*tm.^2))*0.125;
    log_cmb_123 = log_scndtrm_p2 + log_scndtrm_p3 + log_scndtrm_p1;
    scndtrm_p4_var = (exp(6*log_scndtrm_p3)*v_avg.*tm)* g_st;
    log_scndtrm_p4 = scndtrm_p4_var + 0.5 * log(bsxfun(@plus,1,...
        bsxfun(@rdivide,2.026,scndtrm_p4_var)));
    fctr_tmp2 = log_scndtrm_p4 + log_cmb_123 + log_fctr;
    trm_2 = exp(fctr_tmp2- b(1)*tm);
    trm_2(1) = trm_1;
    temp1 = trm_2;
    
    figure(2);
    
    for i=1:100
    plot(log(temp1(10:10000,i)));
    hold on;
    end
    hold off;
    %toc

%end