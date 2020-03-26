function landscape_visualizer
clear;clc;close all;

%load in target data
D = load('data/gill_syn_data_1.mat');

%define visualization colors
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

%define the approximation order to use
%the landscape is a very strong function of i (Taylor order)
%but a weak function of j (Laurent order)
%fig. 2b uses 7,7, but even 7,2 is fine
i=7;
j=7;
F = load(sprintf('landscape/gg_200131_land_scan_1_%i_%i.mat',i,j));

%initialize axes
figure(1);
top = 0.09;
bottom = 0.0;
left = 0.05; right = 0.00;
mid_h = 0.02; mid_v = 0.02;
width = (1-left-right-mid_h)/2;
height = (1-top-bottom-mid_v)/2;

ax_1_pos = [left, bottom+height+mid_v, width, height];
ax_2_pos = [left+width+mid_h, bottom+height+mid_v, width, height];
ax_3_pos = [left, bottom, width, height];
ax_4_pos = [left+width+mid_h, bottom, width, height];

%%%%
% display KL divergence landscape for numerical integration
axes('Position',ax_1_pos);
[~,HH]=contourf(log10(F.KG_),log10(F.B_),...
    reshape(F.data_lik_1,F.n_kg,F.n_b),15);hold on;
set(HH,'LineColor',LINECOL);
scatter(log10(D.kpar(1)), log10(D.S(1)), 40,'r','filled');
colormap(map);
title('Numerical','FontWeight','Normal'); 
ylabel('KL divergence');
h=get(gca);
set(gca,'xtick',[],'ytick',[],'ycolor','none','FontSize',12,'visible','off');
h.YAxis.Label.Color=[0 0 0];
h.YAxis.Label.Visible='on';
h.Title.Color=[0 0 0];
h.Title.Visible='on';


%%%%
% display KL divergence landscape for analytical integration
axes('Position',ax_2_pos);
[~,HH]=contourf(log10(F.KG_),log10(F.B_),...
    reshape(F.data_lik_2,F.n_kg,F.n_b),15);hold on;
set(HH,'LineColor',LINECOL);
scatter(log10(D.kpar(1)), log10(D.S(1)), 40,'r','filled');
colormap(map);
title('Decomposition','FontWeight','Normal'); 
h=get(gca);
set(gca,'xtick',[],'ytick',[],'ycolor','none','FontSize',12,'visible','off')
h.Title.Color=[0 0 0];
h.Title.Visible='on';

%%%%
% display KS distance landscape for numerical integration
axes('Position',ax_3_pos);
[~,HH]=contourf(log10(F.KG_),log10(F.B_),...
    reshape(F.data_emd_1,F.n_kg,F.n_b),15);hold on;
set(HH,'LineColor',LINECOL);
scatter(log10(D.kpar(1)), log10(D.S(1)), 40,'r','filled');
colormap(map);
ylabel('K-S Distance');
h=get(gca);
set(gca,'xtick',[],'ytick',[],'ycolor','none','FontSize',12,'visible','off')
h.YAxis.Label.Color=[0 0 0];
h.YAxis.Label.Visible='on';

%%%%
% display KS distance landscape for analytical integration
axes('Position',ax_4_pos);
[~,HH]=contourf(log10(F.KG_),log10(F.B_),...
    reshape(F.data_emd_2,F.n_kg,F.n_b),15);hold on;
set(HH,'LineColor',LINECOL);
scatter(log10(D.kpar(1)), log10(D.S(1)), 40,'r','filled');
colormap(map);
set(gca,'xtick',[],'ytick',[],'ycolor','none','FontSize',12,'visible','off')

set(gcf,'Position',[488.2000 510.6000 560 251.2000]);
return
