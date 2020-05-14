%profile on;
rw = dir("/home/user/env_EGEL/husn_events/2010-01-01*");
mw = struct2cell(rw);
%length(mw)
%%{
for j=1:length(mw)
    chek = sprintf("%s/%s",mw{(j-1)*7+2},mw{(j-1)*7+1});
    chek
    df_1 = sprintf("%s/%s",chek,'outful_N*.txt');
    alpha = dir(df_1);
    beta = struct2cell(alpha);
    %length(alpha)
    if (~isempty(alpha))
      for p = 1:length(alpha)
        df_3 = sprintf("%s/%s",chek,beta{(p-1)*7+1});
        %df_3
        A = load(df_3);
        b = bartlett(101);
        D = conv(A(1:end-1),b/sum(b),'same');
        pkg load signal
        [q,loc] = findpeaks(D(1001:end));
        E = 1./D;
        [qd,locd] = findpeaks(E(1001:end));
        loc = loc + 1000;
        locd = locd + 1000;
        %length(loc)
        %length(locd)
        x = 1:1:length(D);
        %%plot(x,D); hold on; plot(x(loc),q,'ro*'); plot(x(locd),D(locd),'o'); hold off;
        jmp = zeros(length(loc),1);
        idx = zeros(length(loc),1);
        for i=1:length(loc)
    	    sw = find(locd<loc(i));
	        jmp(i) = 0;
          idx(i) = 0;
	        if (length(sw) > 0)
		        idx(i) = locd(sw(end));
		        jmp(i) = D(loc(i))/D(idx(i));
	        end
        end
        rw = find(jmp > 4);
        cut_p = length(D) - 1;
        if (length(rw) > 0)
	        cut_p = idx(rw(1));
        end

        #fid = fopen(df_3,'w');
        #fprintf(fid,"%e\n",D(1:cut_p)');
        #fprintf(fid,"%f",A(end));
        #fclose(fid);
      end
    end
%profile off;
%T = profile("info");
%profshow(T);
endfor
%}