clc;
close all;
clear

dr_lst = [%dir("/home/user/env_EGEL/husn_events/2010*"); dir("/home/user/env_EGEL/husn_events/2011*");...
    %dir("/home/user/env_EGEL/husn_events/2010-09*"); dir("/home/user/env_EGEL/husn_events/2010-10*");
    %dir("/home/user/env_EGEL/husn_events/2010-11*");... %dir("/home/user/env_EGEL/husn_events/2010-12*");
    %dir("/home/user/env_EGEL/husn_events/2010*"); dir("/home/user/env_EGEL/husn_events/2011*");...
    %dir("/home/user/env_EGEL/husn_events/2012*"); dir("/home/user/env_EGEL/husn_events/2013*"); ...
    %dir("/home/user/env_EGEL/husn_events/2014*"); dir("/home/user/env_EGEL/husn_events/2015*"); ...
    %dir("/home/user/env_EGEL/husn_events/2016*"); ...
    %dir("/home/user/env_EGEL/husn_events/2017*"); ... 
    %dir("/home/user/env_EGEL/husn_events/2018*"); ...
    %dir("/home/user/env_EGEL/husn_events/2019*"); dir("/home/user/env_EGEL/EGEL_events/200*");...
    %dir("/home/user/env_EGEL/Cycnet_events/0*");
    dir("/home/user/env_EGEL/kef_events/2*");];% dir("/home/user/env_EGEL/cornoth_events/2*");];

flag = 0;


for i=1:length(dr_lst)
    i
    chek = sprintf("%s/%s",dr_lst(i).folder,dr_lst(i).name);
    %dr_lst(i).name
    % Accessing the envelope files within the event directory
    df_1 = sprintf("%s/%s",chek,'outful*envwin_*.txt');
    aw = dir(df_1);
    if (~isempty(aw))
        n_fls = length(aw);
        if (n_fls > 0)
            %chek
            for j = 1:n_fls
                df_3 = sprintf("%s/%s",chek,aw(j).name);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fid_i = fopen(df_3,'r');
                temp = textscan(fid_i,'%f');
                fclose(fid_i);
                
                a = split(aw(j).name,'.');
                df_4 = sprintf("%s.bin",a{1});
                df_bin = sprintf("%s/%s",chek,df_4);
                
                fid_b = fopen(df_bin,'w');
                fwrite(fid_b,temp{1},'double');
                fclose(fid_b);
                
                fid_b = fopen(df_bin,'r'); 
                Q = fread(fid_b,'*double'); 
                fclose(fid_b);
                
                diff = Q - temp{1};
                
                if (diff~=0)
                    flag = 1;
                    disp('breaking');
                    break
                end
                delete(df_3);
            end
        end
    end
    if (flag==1)
        disp('error')
        break;
    end
end