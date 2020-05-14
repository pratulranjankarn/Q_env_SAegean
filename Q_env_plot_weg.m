
%%{
clc;
close all;
clearvars;

wrk_path = '/home/user/env_EGEL/env_weg_test';

mkdir(wrk_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 - 2 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fid = fopen('comb_12giwr_weg.txt','r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

sg_lim12 = 5;

freq = 1.5;

n_arr = [A{4} A{5} A{6} A{8} A{9}];


st1 = sprintf('%s/envnodes15g_1-2.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_1-2.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');

factr=1;



un_12 = unique(n_arr(:,1:2),'rows','stable');
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
%ido_ttl = arrayfun(@(x,y)find(n_arr(:,1)==x & n_arr(:,2)==y),un_12(:,1),un_12(:,2),'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = find(n_arr(:,1)==un_12(j,1) & n_arr(:,2)==un_12(j,2));%uint32(ido_ttl{j});
    length(ido)
    if(~isempty(ido))
        v_trk(j) = n_arr(ido(1),3)/factr;
        %length(ido)
        tempg = single(n_arr(ido,4))/factr;
        tempb = single(n_arr(ido,5))/factr;
        s_g = length(tempg);
        if (s_g > sg_lim12)
            mad_g = mad(tempg,1)/0.67;
            mean_g = mean(tempg);
            mad_b = mad(tempb,1)/0.67;
            mean_b = mean(tempb);
            %%{
            for k = 1:100
                %k
                diffg = abs(tempg - mean_g);
                qw = find(diffg < 1.345*mad_g);
                sw = find(diffg >= 1.345*mad_g);
                tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
                mean_g = mean(tempg);
                if (length(qw)==s_g)
                    break;
                end
            end
            
            for k = 1:100
                %k
                diffb = abs(tempb - mean_b);
                qw = find(diffb < 1.345*mad_b);
                sw = find(diffb >= 1.345*mad_b);
                tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
                mean_b = mean(tempb);
                if (length(qw)==s_g)
                    break;
                end
            end
            %%}
            g_f(j) = mean_g;
            b_f(j) = mean_b;
        else
            g_f(j) = NaN;
            b_f(j) = NaN;
        end
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [un_12(:,2) un_12(:,1) log10(Qsc_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid3,'%f %f %f\n',un_arr(q,:)');
fclose(fid3);
un_arr = [un_12(:,2) un_12(:,1) log10(Qi_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid4,'%f %f %f\n',un_arr(q,:)');
fclose(fid4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 - 4 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fid = fopen('comb_24giwr_weg.txt','r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

sg_lim12 = 5;

freq = 1.5;

n_arr = [A{4} A{5} A{6} A{8} A{9}];


st1 = sprintf('%s/envnodes15g_2-4.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_2-4.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');

factr=1;



un_12 = unique(n_arr(:,1:2),'rows','stable');
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
%ido_ttl = arrayfun(@(x,y)find(n_arr(:,1)==x & n_arr(:,2)==y),un_12(:,1),un_12(:,2),'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = find(n_arr(:,1)==un_12(j,1) & n_arr(:,2)==un_12(j,2));%uint32(ido_ttl{j});
    length(ido)
    if(~isempty(ido))
        v_trk(j) = n_arr(ido(1),3)/factr;
        %length(ido)
        tempg = single(n_arr(ido,4))/factr;
        tempb = single(n_arr(ido,5))/factr;
        s_g = length(tempg);
        if (s_g > sg_lim12)
            mad_g = mad(tempg,1)/0.67;
            mean_g = mean(tempg);
            mad_b = mad(tempb,1)/0.67;
            mean_b = mean(tempb);
            %%{
            for k = 1:100
                %k
                diffg = abs(tempg - mean_g);
                qw = find(diffg < 1.345*mad_g);
                sw = find(diffg >= 1.345*mad_g);
                tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
                mean_g = mean(tempg);
                if (length(qw)==s_g)
                    break;
                end
            end
            
            for k = 1:100
                %k
                diffb = abs(tempb - mean_b);
                qw = find(diffb < 1.345*mad_b);
                sw = find(diffb >= 1.345*mad_b);
                tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
                mean_b = mean(tempb);
                if (length(qw)==s_g)
                    break;
                end
            end
            %%}
            g_f(j) = mean_g;
            b_f(j) = mean_b;
        else
            g_f(j) = NaN;
            b_f(j) = NaN;
        end
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [un_12(:,2) un_12(:,1) log10(Qsc_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid3,'%f %f %f\n',un_arr(q,:)');
fclose(fid3);
un_arr = [un_12(:,2) un_12(:,1) log10(Qi_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid4,'%f %f %f\n',un_arr(q,:)');
fclose(fid4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 - 8 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fid = fopen('comb_48giwr_weg.txt','r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

sg_lim12 = 5;

freq = 1.5;

n_arr = [A{4} A{5} A{6} A{8} A{9}];


st1 = sprintf('%s/envnodes15g_4-8.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_4-8.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');

factr=1;



un_12 = unique(n_arr(:,1:2),'rows','stable');
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
%ido_ttl = arrayfun(@(x,y)find(n_arr(:,1)==x & n_arr(:,2)==y),un_12(:,1),un_12(:,2),'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = find(n_arr(:,1)==un_12(j,1) & n_arr(:,2)==un_12(j,2));%uint32(ido_ttl{j});
    length(ido)
    if(~isempty(ido))
        v_trk(j) = n_arr(ido(1),3)/factr;
        %length(ido)
        tempg = single(n_arr(ido,4))/factr;
        tempb = single(n_arr(ido,5))/factr;
        s_g = length(tempg);
        if (s_g > sg_lim12)
            mad_g = mad(tempg,1)/0.67;
            mean_g = mean(tempg);
            mad_b = mad(tempb,1)/0.67;
            mean_b = mean(tempb);
            %%{
            for k = 1:100
                %k
                diffg = abs(tempg - mean_g);
                qw = find(diffg < 1.345*mad_g);
                sw = find(diffg >= 1.345*mad_g);
                tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
                mean_g = mean(tempg);
                if (length(qw)==s_g)
                    break;
                end
            end
            
            for k = 1:100
                %k
                diffb = abs(tempb - mean_b);
                qw = find(diffb < 1.345*mad_b);
                sw = find(diffb >= 1.345*mad_b);
                tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
                mean_b = mean(tempb);
                if (length(qw)==s_g)
                    break;
                end
            end
            %%}
            g_f(j) = mean_g;
            b_f(j) = mean_b;
        else
            g_f(j) = NaN;
            b_f(j) = NaN;
        end
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [un_12(:,2) un_12(:,1) log10(Qsc_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid3,'%f %f %f\n',un_arr(q,:)');
fclose(fid3);
un_arr = [un_12(:,2) un_12(:,1) log10(Qi_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid4,'%f %f %f\n',un_arr(q,:)');
fclose(fid4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8 - 16 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fid = fopen('comb_816giwr_weg.txt','r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

sg_lim12 = 5;

freq = 1.5;

n_arr = [A{4} A{5} A{6} A{8} A{9}];


st1 = sprintf('%s/envnodes15g_8-16.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_8-16.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');

factr=1;



un_12 = unique(n_arr(:,1:2),'rows','stable');
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
%ido_ttl = arrayfun(@(x,y)find(n_arr(:,1)==x & n_arr(:,2)==y),un_12(:,1),un_12(:,2),'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = find(n_arr(:,1)==un_12(j,1) & n_arr(:,2)==un_12(j,2));%uint32(ido_ttl{j});
    length(ido)
    if(~isempty(ido))
        v_trk(j) = n_arr(ido(1),3)/factr;
        %length(ido)
        tempg = single(n_arr(ido,4))/factr;
        tempb = single(n_arr(ido,5))/factr;
        s_g = length(tempg);
        if (s_g > sg_lim12)
            mad_g = mad(tempg,1)/0.67;
            mean_g = mean(tempg);
            mad_b = mad(tempb,1)/0.67;
            mean_b = mean(tempb);
            %%{
            for k = 1:100
                %k
                diffg = abs(tempg - mean_g);
                qw = find(diffg < 1.345*mad_g);
                sw = find(diffg >= 1.345*mad_g);
                tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
                mean_g = mean(tempg);
                if (length(qw)==s_g)
                    break;
                end
            end
            
            for k = 1:100
                %k
                diffb = abs(tempb - mean_b);
                qw = find(diffb < 1.345*mad_b);
                sw = find(diffb >= 1.345*mad_b);
                tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
                mean_b = mean(tempb);
                if (length(qw)==s_g)
                    break;
                end
            end
            %%}
            g_f(j) = mean_g;
            b_f(j) = mean_b;
        else
            g_f(j) = NaN;
            b_f(j) = NaN;
        end
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [un_12(:,2) un_12(:,1) log10(Qsc_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid3,'%f %f %f\n',un_arr(q,:)');
fclose(fid3);
un_arr = [un_12(:,2) un_12(:,1) log10(Qi_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid4,'%f %f %f\n',un_arr(q,:)');
fclose(fid4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16 - 32 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('comb_1632giwr_weg.txt','r');
A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

sg_lim12 = 5;

freq = 1.5;

n_arr = [A{4} A{5} A{6} A{8} A{9}];


st1 = sprintf('%s/envnodes15g_16-32.txt',wrk_path);
st2 = sprintf('%s/envnodes15b_16-32.txt',wrk_path);
fid3 = fopen(st1,'w');
fid4 = fopen(st2,'w');

factr=1;



un_12 = unique(n_arr(:,1:2),'rows','stable');
g_f = zeros(length(un_12),1,'single');
b_f = zeros(length(un_12),1,'single');
v_trk = zeros(length(un_12),1,'single');
%ido_ttl = arrayfun(@(x,y)find(n_arr(:,1)==x & n_arr(:,2)==y),un_12(:,1),un_12(:,2),'UniformOutput',false);
for j=1:length(un_12)
    j
    ido = find(n_arr(:,1)==un_12(j,1) & n_arr(:,2)==un_12(j,2));%uint32(ido_ttl{j});
    length(ido)
    if(~isempty(ido))
        v_trk(j) = n_arr(ido(1),3)/factr;
        %length(ido)
        tempg = single(n_arr(ido,4))/factr;
        tempb = single(n_arr(ido,5))/factr;
        s_g = length(tempg);
        if (s_g > sg_lim12)
            mad_g = mad(tempg,1)/0.67;
            mean_g = mean(tempg);
            mad_b = mad(tempb,1)/0.67;
            mean_b = mean(tempb);
            %%{
            for k = 1:100
                %k
                diffg = abs(tempg - mean_g);
                qw = find(diffg < 1.345*mad_g);
                sw = find(diffg >= 1.345*mad_g);
                tempg(sw) = tempg(sw) .* 1.345 .* mad_g./diffg(sw);
                mean_g = mean(tempg);
                if (length(qw)==s_g)
                    break;
                end
            end
            
            for k = 1:100
                %k
                diffb = abs(tempb - mean_b);
                qw = find(diffb < 1.345*mad_b);
                sw = find(diffb >= 1.345*mad_b);
                tempb(sw) = tempb(sw) .* 1.345 .* mad_b./diffb(sw);
                mean_b = mean(tempb);
                if (length(qw)==s_g)
                    break;
                end
            end
            %%}
            g_f(j) = mean_g;
            b_f(j) = mean_b;
        else
            g_f(j) = NaN;
            b_f(j) = NaN;
        end
    end
end

Qsc_inv = g_f.*v_trk/(2*pi*freq);
Qi_inv = b_f/(2*pi*freq);

un_arr = [un_12(:,2) un_12(:,1) log10(Qsc_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid3,'%f %f %f\n',un_arr(q,:)');
fclose(fid3);
un_arr = [un_12(:,2) un_12(:,1) log10(Qi_inv)];
q = ~isnan(un_arr(:,3));
fprintf(fid4,'%f %f %f\n',un_arr(q,:)');
fclose(fid4);

clear n_arr un_12 g_f b_f v_trk un_arr tempg tempb diffg diffb qw sw