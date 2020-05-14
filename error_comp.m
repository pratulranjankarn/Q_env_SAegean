clc;
close all;
clearvars;

wkdir = '/home/user/env_EGEL/env_gauss_7.2,3.6,1.8,.90_pkdbl';


abs_pdir = sprintf('%s/env_err_1-2.txt',wkdir);
fid = fopen(abs_pdir,'r');
A = textscan(fid,'%s %f %f %f');
fclose(fid);
% fid = fopen('envweg_err_1-2.txt','r');
% B = textscan(fid,'%s %f %f');
% fclose(fid);

error_prom  = A{3};
error_prom(error_prom > 1) = [];
error_prom(isnan(error_prom)) = [];

error_sum  = A{4};
error_sum(error_sum > 100) = [];
error_sum(isnan(error_sum)) = [];

%n_a = round(length(A{2})/100);

t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_prom));
t1_st_f = sprintf("%s = %f",t1_st,std(error_prom));
% t2_f = sprintf("%s = %f",t2,mean(B{2}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{2}));

FigH = figure('Position', get(0, 'Screensize')); ...
subplot(4,2,1);
%histogram(error_prom,'BinWidth',0.01,'FaceColor',[0 0.5 0.5]);
h = histfit(error_prom,100,'exponential'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{2}));
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [0.85 0.5];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 1]);
    text(x_pos(1),y_pos(2)*2.5/5,'1-2 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ... 
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('S-window % error'); hold off;
    
    pd = fitdist(error_prom,'exponential')

t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_sum));
t1_st_f = sprintf("%s = %f",t1_st,std(error_sum));
% t2_f = sprintf("%s = %f",t2,mean(B{3}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{3}));

%x_pos = max(max(log10(error_prom)),max(log10(B{3}))) * 3/4;

n_ft = length(min(error_sum):1:max(error_sum));

subplot(4,2,2);
%histogram(error_sum,'BinWidth',1,'FaceColor',[0 0.5 0.5]); hold on;...
h = histfit(error_sum,n_ft,'normal'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{3})); ...
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [85 55];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 100]);
    text(x_pos(1),y_pos(2)*2.5/5,'1-2 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('Total % error'); hold off;
% fig = gcf;
% fig.PaperPositionMode = 'manual';
% fig.PaperPosition = [.25 .25 8 6];
% print('env_comp_1-2 Hz','-dpdf','-r0')
%%

abs_pdir = sprintf('%s/env_err_2-4.txt',wkdir);
fid = fopen(abs_pdir,'r');
A = textscan(fid,'%s %f %f %f');
fclose(fid);
% fid = fopen('envweg_err_2-4.txt','r');
% B = textscan(fid,'%s %f %f');
% fclose(fid);

error_prom  = A{3};
error_prom(error_prom > 1) = [];
error_prom(isnan(error_prom)) = [];

error_sum  = A{4};
error_sum(error_sum > 100) = [];
error_sum(isnan(error_sum)) = [];

%n_a = round(length(A{2})/100);

t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_prom));
t1_st_f = sprintf("%s = %f",t1_st,std(error_prom));
% t2_f = sprintf("%s = %f",t2,mean(B{2}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{2}));

subplot(4,2,3);
%histogram(error_prom,'BinWidth',0.01,'FaceColor',[0 0.5 0.5]);
h = histfit(error_prom,100,'exponential'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{2}));
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [0.85 0.5];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 1]);
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ... 
    text(x_pos(1),y_pos(2)*2.5/5,'2-4 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('S-window % error'); hold off;
    
    pd = fitdist(error_prom,'exponential')


t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_sum));
t1_st_f = sprintf("%s = %f",t1_st,std(error_sum));
% t2_f = sprintf("%s = %f",t2,mean(B{3}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{3}));

%x_pos = max(max(log10(error_prom)),max(log10(B{3}))) * 3/4;

n_ft = length(min(error_sum):1:max(error_sum));

subplot(4,2,4);
%histogram(error_sum,'BinWidth',1,'FaceColor',[0 0.5 0.5]); hold on;...
h = histfit(error_sum,n_ft,'normal'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{3})); ...
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [85 55];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 100]);
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(1),y_pos(2)*2.5/5,'2-4 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('Total % error'); hold off;
% fig = gcf;
% fig.PaperPositionMode = 'manual';
% fig.PaperPosition = [.25 .25 8 6];
% print('env_comp_2-4 Hz','-dpdf','-r0')
%%
abs_pdir = sprintf('%s/env_err_4-8.txt',wkdir);
fid = fopen(abs_pdir,'r');
A = textscan(fid,'%s %f %f %f');
fclose(fid);
% fid = fopen('envweg_err_4-8.txt','r');
% B = textscan(fid,'%s %f %f');
% fclose(fid);

error_prom  = A{3};
error_prom(error_prom > 1) = [];
error_prom(isnan(error_prom)) = [];

error_sum  = A{4};
error_sum(error_sum > 100) = [];
error_sum(isnan(error_sum)) = [];

%n_a = round(length(A{2})/100);

t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_prom));
t1_st_f = sprintf("%s = %f",t1_st,std(error_prom));
% t2_f = sprintf("%s = %f",t2,mean(B{2}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{2}));

subplot(4,2,5);
%histogram(error_prom,'BinWidth',0.01,'FaceColor',[0 0.5 0.5]);
h = histfit(error_prom,100,'exponential'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{2}));
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [0.85 0.5];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 1]);
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ... 
    text(x_pos(1),y_pos(2)*2.5/5,'4-8 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('S-window % error'); hold off;
    
    pd = fitdist(error_prom,'exponential')


t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_sum));
t1_st_f = sprintf("%s = %f",t1_st,std(error_sum));
% t2_f = sprintf("%s = %f",t2,mean(B{3}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{3}));

%x_pos = max(max(log10(error_prom)),max(log10(B{3}))) * 3/4;

n_ft = length(min(error_sum):1:max(error_sum));

subplot(4,2,6);
%histogram(error_sum,'BinWidth',1,'FaceColor',[0 0.5 0.5]); hold on;...
h = histfit(error_sum,n_ft,'normal'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts'); ...
%     histogram(log10(B{3})); ...
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [85 55];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 100]);
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(1),y_pos(2)*2.5/5,'4-8 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('Total % error'); hold off;
% fig = gcf;
% fig.PaperPositionMode = 'manual';
% fig.PaperPosition = [.25 .25 8 6];
% print('env_comp_4-8 Hz','-dpdf','-r0')
%%
abs_pdir = sprintf('%s/env_err_8-16.txt',wkdir);
fid = fopen(abs_pdir,'r');
A = textscan(fid,'%s %f %f %f');
fclose(fid);
% fid = fopen('envweg_err_8-16.txt','r');
% B = textscan(fid,'%s %f %f');
% fclose(fid);

error_prom  = A{3};
error_prom(error_prom > 1) = [];
error_prom(isnan(error_prom)) = [];

error_sum  = A{4};
error_sum(error_sum > 100) = [];
error_sum(isnan(error_sum)) = [];

%n_a = round(length(A{2})/100);

t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_prom));
t1_st_f = sprintf("%s = %f",t1_st,std(error_prom));
% t2_f = sprintf("%s = %f",t2,mean(B{2}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{2}));

% pd = fitdist(error_prom,'exponential');
% x_values = -1:0.01:1;
% pd2 = pdf(pd,x_values);

subplot(4,2,7);
%histogram(error_prom,'BinWidth',0.01,'FaceColor',[0 0.5 0.5]);
h = histfit(error_prom,100,'exponential'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{2}));
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [0.85 0.5];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 1]);
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(1),y_pos(2)*2.5/5,'8-16 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('S-window % error'); hold off;
    
    pd = fitdist(error_prom,'exponential')

t1 = '\mu\fontsize{14}';
t1_st = '\sigma\fontsize{14}';
% t2 = '\mu_{WinAvg.}\fontsize{14}';
% t2_st = '\sigma_{WinAvg.}\fontsize{14}';

t1_f = sprintf("%s = %f",t1,mean(error_sum));
t1_st_f = sprintf("%s = %f",t1_st,std(error_sum));
% t2_f = sprintf("%s = %f",t2,mean(B{3}));
% t2_st_f = sprintf("%s = %f",t2_st,std(B{3}));

%x_pos = max(max(log10(error_prom)),max(log10(B{3}))) * 3/4;

n_ft = length(min(error_sum):1:max(error_sum));

subplot(4,2,8); 
%histogram(error_sum,'BinWidth',1,'FaceColor',[0 0.5 0.5]); hold on;...
h = histfit(error_sum,n_ft,'normal'); hold on; ...
    h(1).FaceColor = [0 0.5 0.5];
%     ylabel('Counts');...
%     histogram(log10(B{3})); ...
%     xlabel('% Error'); ...
    set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');...
    %x_pos = get(gca,'XLim')*0.6
    x_pos = [85 55];
    y_pos = get(gca,'YLim')
    set(gca,'XLim',[0 100]);
    text(x_pos(2),y_pos(2)*2.5/5,t1_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(2),y_pos(2)*1.5/5,t1_st_f,'FontSize',14,'FontWeight','bold'); ...
    text(x_pos(1),y_pos(2)*2.5/5,'8-16 Hz','Color','red','FontSize',14,'FontWeight','bold'); ...
%     text(x_pos(2),y_pos(2)*3/5,t2_f,'FontSize',12); ...
%     text(x_pos(2),y_pos(2)*2.5/5,t2_st_f,'FontSize',12); ...
    legend('Total % error'); hold off;
% fig = gcf;
% fig.PaperPositionMode = 'manual';
% fig.PaperPosition = [.25 .25 8 6];
% print('env_comp_8-16 Hz','-dpdf','-r0')

%%

h4 = suplabel('Counts','y');
set(h4,'Position',[0.110 0.1000 0.7750 0.8150],'FontSize',20,'FontWeight','bold');

h5 = suplabel('% Error','x');
set(h5,'Position',[0.120 0.1000 0.7750 0.8150],'FontSize',18,'FontWeight','bold');

F    = getframe(FigH);
imwrite(F.cdata, 'env_err_comp_pkdbl', 'png')