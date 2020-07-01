function supp_fig_3
clear;clc;close all;
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
BATCH = 12;
SCAN = 24;

analyt_col = [217 116 55]/255;
LINECOL = 'none';

p1=[227,172,82]/255;
p2=[252,222,164]/255;
p3=[90,180,172]/255;
n=8;
map = NaN(n*2,3);
for i = 1:3
    map(1:n,i) = linspace(p1(i),p2(i),n);
    map((n+1):end,i) = linspace(p2(i),p3(i),n);
end

c1 = 226/266;
c2 = 177/266;
n=8;
map_grey = repmat(linspace(c2,c1,n)',[1,3]);
figure(1);
N_horz = 6; N_vert = 4;

DAT_name_1 = 'KL divergence';
DAT_name_2 = 'chf divergence';

fig1=figure(1);
ha1=tight_subplot(N_vert,N_horz);
fig1.Position=[193.8000 184.2000 911.2000 522.4000];

fig2=figure(2);
ha2=tight_subplot(N_vert,N_horz);
fig2.Position=[193.8000 184.2000 911.2000 522.4000];
for DATANUM = 1:N_horz*N_vert
    D = load(sprintf('data/gg_200618_sim_%i_%i.mat',BATCH,DATANUM));
    F = load(sprintf('landscape/gg_200625_land_scan_%i_%i_%i_%i.mat',SCAN,DATANUM,7,1));
    
%     DATA_TO_PLOT = F.data_lik_1;
    DATA_TO_PLOT = F.data_lik_2;
%     DATA_TO_PLOT = log10(F.data_chf_1);
%     DATA_TO_PLOT = log10(F.data_chf_2); %this is analytically evaluated chf diff
    DATA_THRESHOLD = quantile(DATA_TO_PLOT,0.05);
%     
    axes(ha1(DATANUM));
    [~,HH]=contourf(log10(F.KG_),log10(F.B_),...
        reshape(DATA_TO_PLOT,F.n_kg,F.n_b),100);hold on;
    set(HH,'LineColor',LINECOL);
    scatter(log10(D.kpar(1)), log10(D.S(1)), 40,'r','filled');
    colormap(map);

    h=get(gca);
    set(gca,'xtick',[],'ytick',[],...
        'ycolor','none','FontSize',12,'visible','off','box','on');
    h.YAxis.Label.Color=[0 0 0];
    h.YAxis.Label.Visible='on';
    h.Title.Color=[0 0 0];
    h.Title.Visible='on';

    axes(ha2(DATANUM));
    [~,HH]=contourf(log10(F.KG_),log10(F.B_),...
        reshape(DATA_TO_PLOT>DATA_THRESHOLD,F.n_kg,F.n_b),100);hold on;
    set(HH,'LineColor',LINECOL);
    scatter(log10(D.kpar(1)), log10(D.S(1)), 40,'r','filled');
    colormap(map_grey);

    h=get(gca);
    set(gca,'xtick',[],'ytick',[],...
        'ycolor','none','FontSize',12,'visible','off');
    h.YAxis.Label.Color=[0 0 0];
    h.YAxis.Label.Visible='on';
    h.Title.Color=[0 0 0];
    h.Title.Visible='on';
end
return
