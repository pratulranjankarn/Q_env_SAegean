clc;
close all;
clearvars;


fid = fopen('env_WR_1-2.txt','r');
A1 = textscan(fid,"%f %f %f %f");
fclose(fid);

fid = fopen('env_WR_2-4.txt','r');
A2 = textscan(fid,"%f %f %f %f");
fclose(fid);

fid = fopen('env_WR_4-8.txt','r');
A3 = textscan(fid,"%f %f %f %f");
fclose(fid);

fid = fopen('env_WR_8-16.txt','r');
A4 = textscan(fid,"%f %f %f %f");
fclose(fid);

fid = fopen('env_WR_16-32.txt','r');
A5 = textscan(fid,"%f %f %f %f");
fclose(fid);

A = [cell2mat(A1); cell2mat(A2); cell2mat(A3); cell2mat(A4); cell2mat(A5)];

evnttl = max(unique(A(:,1)));

sttl = 187*5;

evn_idx = A(:,1);
st_idx = A(:,2);
freq = A(:,3);

dta = A(:,4);
evn_idx(isinf(log(dta))) = [];
st_idx(isinf(log(dta))) = [];
freq(isinf(log(dta))) = [];
dta(isinf(log(dta))) = [];

n_a = length(st_idx);

w_R = 1;

parm1 = logspace(3,6,10);
parm2 = linspace(1,1,10);

modl = zeros(n_a,1);

post_prob = zeros(1,100);
param_set = zeros(100,2);
m=1;

for i = 1:10
    for j = 1:10
        param_set(m,:) = [parm1(i) parm2(j)];
        vrs1 = ones(evnttl,1)*parm1(i);
        vrs2 = ones(sttl,1)*parm2(j);
        for k=1:n_a
            modl(k) = log(vrs1(evn_idx(k))) + log(vrs2(st_idx(k)));
        end
        post_prob(m) = -sum((log(dta)-modl).^2);% - w_R*(1-(prod(vrs2).^(1/sttl))).^2;
        m=m+1;
    end
end

m = find(post_prob==max(post_prob))



%%{

vrs1 = ones(evnttl,1).*rand(evnttl,1)*param_set(m(1),1);
vrs2 = ones(sttl,1).*rand(sttl,1)*param_set(m(1),2);

n_iter = 500000;
count=0;
n_burn = 300000;

for j=1:n_a
        modl(j) = log(vrs1(evn_idx(j))) + log(vrs2(st_idx(j)));
end
post_prob = -sum((log(dta)-modl).^2) - w_R*(1-(prod(vrs2).^(1/sttl))).^2;

vrs1_f = 0;
vrs2_f = 0;


for i = 1:n_iter
    [vrs1,vrs2,count,post_prob,modl] = lin_fn(vrs1,vrs2,dta,evn_idx,st_idx,modl,sttl,evnttl,freq,post_prob,i,n_burn,count,w_R);
    %count
    post_prob
    i
    if (count > 1)
       vrs1_f = (vrs1_f * (count-1) + vrs1)/count;
       vrs2_f = (vrs2_f * (count-1) + vrs2)/count;
    end
end
%}

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
    dir("/home/user/env_EGEL/Cycnet_events/0*");];

mw = struct2cell(dr_lst);

W_f = log10(vrs1);
fid = fopen('W_f.txt','w');
for i =1:(length(W_f))
fprintf(fid,'%s %f\n',dr_lst(i).name,W_f(i));
end
fclose(fid);

ind1 = round(sttl/5);
ind2 = round(2*sttl/5);
ind3 = round(3*sttl/5);
ind4 = round(4*sttl/5);
        
site_f = log10(vrs2);
fid = fopen('stlist_f.txt','r');
st_db = textscan(fid,'%s %f %f');
fclose(fid);

fid = fopen('R_f_1-2.txt','w');
R_f = [st_db{3} st_db{2} site_f(1:ind1)];
fprintf(fid,'%f %f %f\n',R_f');
fclose(fid);

fid = fopen('R_f_2-4.txt','w');
R_f = [st_db{3} st_db{2} site_f((ind1+1):ind2)];
fprintf(fid,'%f %f %f\n',R_f');
fclose(fid);

fid = fopen('R_f_4-8.txt','w');
R_f = [st_db{3} st_db{2} site_f((ind2+1):ind3)];
fprintf(fid,'%f %f %f\n',R_f');
fclose(fid);

fid = fopen('R_f_8-16.txt','w');
R_f = [st_db{3} st_db{2} site_f((ind3+1):ind4)];
fprintf(fid,'%f %f %f\n',R_f');
fclose(fid);

fid = fopen('R_f_16-32.txt','w');
R_f = [st_db{3} st_db{2} site_f((ind4+1):end)];
fprintf(fid,'%f %f %f\n',R_f');
fclose(fid);

function [vrs1,vrs2,count,post_prob,modl] = lin_fn(vrs1,vrs2,dta,evn_idx,st_idx,modl,sttl,evnttl,freq,post_prob,iter,n_burn,count,w_R)
        std_vrs1 = 50;
        std_vrs2 = 1;
        vrs1_tmp = vrs1;
        vrs2_tmp = vrs2;
        chus = randi(2);
        if chus == 1
            ind = single(randi(evnttl));
            temp1 = vrs1_tmp(ind);
            temp2 = single(trandn((1000-temp1)/std_vrs1,(10^6-temp1)/std_vrs1));
            vrs1_tmp(ind) = temp1 + temp2 * std_vrs1;
            rw = find(evn_idx==ind);
        else
            ind = single(randi(sttl));
            temp1 = vrs2_tmp(ind);
            temp2 = single(trandn((0.001-temp1)/std_vrs2,(1000-temp1)/std_vrs2));
            vrs2_tmp(ind) = temp1 + temp2 * std_vrs2;
            if ind < round(sttl/5)
                frq_tmp = 1.5;
            elseif (ind > round(sttl/5)) && (ind <= round(2*sttl/5))
                frq_tmp = 3;
            elseif (ind > round(2*sttl/5)) && (ind <= round(3*sttl/5))
                frq_tmp = 6;
            elseif (ind > round(3*sttl/5)) && (ind <= round(4*sttl/5))
                frq_tmp = 12;
            else
                frq_tmp = 24;
            end
            rw = find(st_idx==ind & freq==frq_tmp);
        end
        modl(rw) = log(vrs1_tmp(evn_idx(rw))) + log(vrs2_tmp(st_idx(rw)));
        ind1 = round(sttl/5);
        ind2 = round(2*sttl/5);
        ind3 = round(3*sttl/5);
        ind4 = round(4*sttl/5);
        post_prob1 = -sum((log(dta)-modl).^2) - w_R*(1-(prod(vrs2(1:ind1)).^(1/ind1))).^2 ...
            - w_R*(1-(prod(vrs2((ind1+1):ind2)).^(1/(ind2-ind1)))).^2 ...
            - w_R*(1-(prod(vrs2((ind2+1):end)).^(1/(ind3-ind2)))).^2 ...
            - w_R*(1-(prod(vrs2((ind3+1):end)).^(1/(ind4-ind3)))).^2 ...
            - w_R*(1-(prod(vrs2((ind4+1):end)).^(1/(sttl-ind4)))).^2;
        paccept = post_prob1 - post_prob;
        accep_prob1 = min(1,exp(paccept));
        accep_prob1
        if (accep_prob1 == 1 || accep_prob1 > rand(1))
            vrs1 = vrs1_tmp;
            vrs2 = vrs2_tmp;
            post_prob = post_prob1;
            if iter > n_burn
                count = count+1;
            end
        end
end