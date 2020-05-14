function error_sum = g_func_wrapper2(g_st,w_len,E,tm,v_avg,r,sz,t_end)
%sum_tmp=0;

temp1 = zeros(sz,1,'single');
g_st = single(g_st);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_len = single((w_len)/100);

tcons = single(r/v_avg);


tm_dm = single(linspace(tcons,t_len+tcons,(w_len)+1));
%length(tm_dm)
scndtrm_p5 = zeros(1,w_len+1,'single');

fctr = exp(-v_avg*tm_dm*g_st(1));
frstrm = zeros(1,w_len+1,'single');
[~,qw] = min(tm_dm-r/v_avg,[],2);
frstrm(qw) = 1/(4*pi*r^2);
scndtrm_p1 = exp(-1.5*(log(4*pi*v_avg)-log(3*g_st(1))));
scndtrm_p2 = exp((-1.5)*log(tm_dm));
alpha = r./(v_avg.*tm_dm);
beta = (1 - alpha.^2);
beta(qw) = 0;
scndtrm_p3 = exp(log(beta)*(1/8));
chek = v_avg*tm_dm-r;
scndtrm_p5(chek > 0) = 1;
scndtrm_p5(chek == 0) = 0.5;
scndtrm_p5(chek < 0) = 0;
%scndtrm_p5 = heaviside(v_avg*tm_dm-r);
scndtrm_p4_var = exp(log(v_avg*tm_dm*g_st(1)) + 0.75*log(beta));
scndtrm_p4_var(qw) = 0;
scndtrm_p4 = exp(scndtrm_p4_var + 0.5*log(1 + 2.026./scndtrm_p4_var));
scndtrm_p4(qw) = 0;
%max(fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5))
temp1(1:w_len+1) = fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5);%/((w_len-100)+1);


tm_dm = single(linspace(t_len+0.01+tcons,t_end+tcons,(sz-w_len-1)));
%length(tm_dm)
scndtrm_p5 = zeros(1,(sz-w_len-1),'single');

fctr = exp(-v_avg*tm_dm*g_st(1));
frstrm = single(0);
scndtrm_p1 = (4*pi*v_avg/(3*g_st(1)))^(-1.5);
scndtrm_p2 = exp((-1.5)*log(tm_dm));
scndtrm_p3 = exp(log(1 - r^2./(v_avg^2*tm_dm.^2))*(0.125));
scndtrm_p4_var = v_avg*tm_dm*g_st(1).*(1-r^2./(v_avg^2*tm_dm.^2)).^0.75;
scndtrm_p4 = exp(scndtrm_p4_var).* sqrt(1 + 2.026./scndtrm_p4_var);
chek = v_avg*tm_dm-r;
scndtrm_p5(chek > 0) = 1;
scndtrm_p5(chek == 0) = 0.5;
scndtrm_p5(chek < 0) = 0;
temp1(w_len+2:sz) = fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_f1 = log(temp1);

G_f(1) = mean(G_f1(1:w_len));
G_f(2:(sz-w_len)) = G_f1((w_len+1):end-1);

%size(G_f)
%size(E)
%size(tm)

X = [ones(length(G_f),1,'single') -tm'];

%{

lb = [0 0];
options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
%options = optimoptions(options,'SubproblemAlgorithm','factorization');
%options = optimoptions(options,'FunctionTolerance',1.000000e-10,'OptimalityTolerance',1e-34);
[b,resnorm,residual,exitflag,output] = ...
              lsqlin(double(X),double((E-G_f)'),[],[],[],[],lb,[],[],options);
%}

b = X\(E-G_f);
sum_tmp = ((E-G_f)-X*b).^2;
%E_calc = G_f+b(1,j)-b(2,j).*tm';
%t_fnl = linspace(-1.0,tm(end),sz-1);
%figure(j)
%scatter(tm(1),E(1)); hold on; plot(tm,E); %legend(string(r));
%hold on; scatter(tm(1),E_calc(1)); plot(tm,E_calc); plot(t_fnl,log(E_arr(1:end-1))); hold off;
error_sum = double(sum_tmp);
%{
if (b(2) < 0)
    error_sum = inf;
end
%}
if (b(2)==0) || (isnan(b(2)))
b = [];
end
if (~isreal(b)) || (isempty(b))
    error_sum = inf;
end
%clearvars('-except','error_sum');
end