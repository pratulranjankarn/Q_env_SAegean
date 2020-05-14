function error_sum = g_func_wrapper_weg(g_st,w_len,Earset,tm_set,v_set,r_set,sz_set,n_fls)
g_st = single(g_st);
G_f_set = cell(1,n_fls);

sum_tmp = cell(1,n_fls);

%Cnst = cell(1,n_fls);

b_inv = zeros(n_fls,2,'single');

for i=1:n_fls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_avg = v_set(i);
r = r_set(i);
sz = sz_set(i);

t_end = (sz-100)/100;

temp1 = zeros(1,sz,'single');

t_len = single((w_len-100)/100);

tcons = single(r/v_avg);

tm_dm = single(linspace(tcons,t_len+tcons,(w_len-100)+1));
%length(tm_dm)
scndtrm_p5 = zeros(1,w_len-99,'single');

fctr = exp(-v_avg*tm_dm*g_st(i));
frstrm = zeros(1,length(tm_dm));
[~,qw] = min(tm_dm-r/v_avg,[],2);
frstrm(qw) = 1/(4*pi*r^2);
scndtrm_p1 = exp(-1.5*(log(4*pi*v_avg)-log(3*g_st(i))));
scndtrm_p2 = exp((-1.5)*log(tm_dm));
alpha = r./(v_avg.*tm_dm);
beta = (1 - alpha.^2);
beta(qw) = 0;
scndtrm_p3 = exp(log(beta)*(1/8));
chek = v_avg*tm_dm-r;
scndtrm_p5(chek > 0) = 1;
scndtrm_p5(chek == 0) = 0.5;
scndtrm_p5(chek < 0) = 0;
scndtrm_p4_var = exp(log(v_avg*tm_dm*g_st(i)) + 0.75*log(beta));
scndtrm_p4_var(qw) = 0;
scndtrm_p4 = exp(scndtrm_p4_var + 0.5*log(1 + 2.026./scndtrm_p4_var));
scndtrm_p4(qw) = 0;
temp1(1:(w_len-99)) = fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5);


tm_dm = linspace(t_len+0.01+tcons,t_end+tcons,((sz-w_len)+99));
scndtrm_p5 = zeros(1,(sz-w_len+99),'single');

fctr = exp(-v_avg*tm_dm*g_st(i));
frstrm = 0;
scndtrm_p1 = (4*pi*v_avg/(3*g_st(i)))^(-1.5);
scndtrm_p2 = exp((-1.5)*log(tm_dm));
scndtrm_p3 = exp(log(1 - r^2./(v_avg^2*tm_dm.^2))*(0.125));
scndtrm_p4_var = v_avg*tm_dm*g_st(i).*(1-r^2./(v_avg^2*tm_dm.^2)).^0.75;
scndtrm_p4 = exp(scndtrm_p4_var).* sqrt(1 + 2.026./scndtrm_p4_var);
chek = v_avg*tm_dm-r;
scndtrm_p5(chek > 0) = 1;
scndtrm_p5(chek == 0) = 0.5;
scndtrm_p5(chek < 0) = 0;
temp1((w_len-98):sz) = fctr .* (frstrm + scndtrm_p1 .* scndtrm_p2 .* scndtrm_p3 .* scndtrm_p4 .* scndtrm_p5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_f1 = log(temp1);

G_f_set{i}(1) = mean(G_f1(1:w_len));
G_f_set{i}(2:(sz-w_len)) = G_f1((w_len+1):end-1);

X = [ones((sz-w_len),1,'single') -(tm_set{i})'];


%{
lb = [0 1e-6];
%ub = [10 0.5];
options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
options.Display = 'off';
options = optimoptions(options,'SubproblemAlgorithm','factorization');
%options = optimoptions(options,'FunctionTolerance',1.000000e-10,'OptimalityTolerance',1e-34);
[b,resnorm,residual,exitflag,output] = ...
    lsqlin(double(X),double((Earset{i}-G_f_set{i})'),[],[],[],[],lb,[],[],options);
%}

b = X\(Earset{i}-G_f_set{i})';

b_inv(i,:) = b;

%size(Earset{i})
%size(G_f_set{i})
%size(tm_set{i})
%size(b)

sum_tmp{i} = (((Earset{i}-G_f_set{i})'-X*b).^2)';

%Cnst{i} = ones(1,(sz-w_len),'single') * R(i);
end

%{

G_f = cell2mat(G_f_set);

Cnst_ar = cell2mat(Cnst);


X = [ones(length(G_f),1,'single') -tm_set'];

lb = [0 0];
%ub = [10 0.5];
options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
options.Display = 'off';
%options = optimoptions(options,'SubproblemAlgorithm','factorization');
%options = optimoptions(options,'FunctionTolerance',1.000000e-10,'OptimalityTolerance',1e-34);
[b,resnorm,residual,exitflag,output] = ...
    lsqlin(double(X),double((Earset-G_f)'),[],[],[],[],lb,[],[],options);

            
sum_tmp = ((Earset-G_f)'-X*b).^2;
%E_calc = G_f+b(1)-b(2).*tm';
%t_fnl = linspace(-1.0,tm(end),sz-1);
%figure(j)
%scatter(tm(1),E(1)); hold on; plot(tm,E); %legend(string(r));
%hold on; scatter(tm(1),E_calc(1)); plot(tm,E_calc); plot(t_fnl,log(E_arr(1:end-1))); hold off;
%b
%}
error_sum = double(cell2mat(sum_tmp));
if (~isreal(b)) || (isempty(b)) %|| (~isreal(g_st))
    error_sum = inf;
end
%clearvars('-except','error_sum');
end