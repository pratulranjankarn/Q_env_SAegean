function error_sum = g_func_wrpr_mtrx(g_st,w_len,E,tm,v_avg,r,sz)

% g_st = single(logspace(-5,-1,100));
% v_avg = single(3.3373);
% r = single(40.830);
% sz = single(10000);
% w_len = single(1000);
% t_end = single(100);

    g_st = single(g_st);
    
    n_g = length(g_st);

    %G_f = ones((sz-w_len),length(g_st),'single');
    
    t_len = single(sz/100);

    tcons = single(r/v_avg);
    
    
    tm_dm = (tcons:0.01:(t_len+tcons-0.01))';
    
    log_fctr = -v_avg*tm_dm*g_st;
    frstrm = 1/(4*pi*r^2);
    fctr_tmp1 = log_fctr(1,:) + log(frstrm);
    lg_trm_1 = fctr_tmp1;
    
    log_scndtrm_p1 = -1.5*(log(4*pi*v_avg)-log(3*g_st));
    log_scndtrm_p2 = (-1.5)*log(tm_dm);
    log_scndtrm_p3 = log(1 - r^2./(v_avg^2*tm_dm.^2))*0.125;
    
    log_cmb_123 = log_scndtrm_p2 + log_scndtrm_p3 + log_scndtrm_p1;
    scndtrm_p4_var = (exp(6*log_scndtrm_p3)*v_avg.*tm_dm)* g_st;
    log_scndtrm_p4 = scndtrm_p4_var + 0.5 * log(bsxfun(@plus,1,...
        bsxfun(@rdivide,2.026,scndtrm_p4_var)));
    fctr_tmp2 = log_scndtrm_p4 + log_cmb_123 + log_fctr;
    lg_trm_2 = fctr_tmp2;
    lg_trm_2(1,1:n_g) = lg_trm_1;
    lg_temp1 = lg_trm_2;
    
    %length(temp1(isnan(temp1)))
    %length(temp1(isinf(temp1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lG_f1 = bsxfun(@minus,lg_temp1,20.723266);

alpha1 = log(mean(exp(lG_f1(1:w_len+1,1:n_g))));
alpha2 = lG_f1((w_len+2):sz,1:n_g); G_f = cat(1,alpha1,alpha2);

X = [ones((sz-w_len),1,'single') -tm'];

% size(G_f1)
% size(G_f)
% size(X)
% size(E)
% size(E-G_f)

b = X\(E-G_f);

%{

b = zeros(2,100);

lb = [0 1e-3];
ub = [inf 1e-1];

for i=1:100
b(:,i) = lsqlin(double(X),double((E-G_f(:,1))),[],...
    [],[],[],lb,ub,[],options);
end
%}

idx = b(2,:) < 0;
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
error_sum = double((sum(sum_tmp(2:end,:)) + sum_tmp(1,:)*w_len));

error_sum(idx) = inf;
%%{
if (~isreal(b)) || (isempty(b)) %|| (any(b(2,:)<0))
    error_sum = inf;
end
%}
%clearvars('-except','error_sum');
end