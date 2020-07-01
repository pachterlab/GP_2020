function fig_2a
clear;clc;close all;

%a e s t h e t i c s
line_width = 0.8;
line_color = 'r';
face_color = 0.4*[1 1 1];
edge_color = 'none';

%initialize physiology
beta = 1.1;
gamma = beta;
kinit= 2;
bs=10;

%initialize simulation settings
nCells = 1e5;
rng(floor(log(4815162342)));
S = [bs 0; -1 1; 0 -1];
nT = 4;
kpar = [kinit,beta,gamma];
Time_max = 1.2;
tvec = linspace(0,Time_max,nT);
t_matrix = repmat(tvec,nCells,1);

tic
X = gg_200110_gillespie_geom_1(kpar,t_matrix,S,nCells);
toc

%initialize range of approximation settings
Tmax = 8;
Lmax = 4;


f=figure;
set(f,'Position',[636 200 604 230]);

for NN = 1:(Tmax*Lmax)
    N_approx_taylor = mod(NN-1,Tmax)+1;
    N_approx_laurent = floor((NN-1)/Tmax)+1;

    tvec(end) = 1.2;
    
    %plot simulation results (mature marginal)
    figure(1)
    subplot(Lmax,Tmax,NN)
    histogram(X(:,end,2),'BinMethod','integers','Normalization','pdf',...
        'FaceColor',face_color,'EdgeColor',edge_color); hold on;
    set(gca,'xtick',[],'ytick',[],'box','off');

    M = max(X(:,end,1))+1;
    N = max(X(:,end,2))+1;

    %compute approximation results
    Pa_marg =  gg_200325_analyt_geom_tdep_vec_31(kinit,bs,gamma,M,N,tvec(end),...
        'mature',N_approx_taylor,N_approx_laurent);
    
    %plot approximation results
    plot(0:(N-1),Pa_marg,strcat(line_color,'-'),'LineWidth',line_width);
    
    %if desired, output approximation order as title of each plot
%     title(sprintf('T=%i, L=%i',N_approx_taylor,N_approx_laurent),...
%         'FontWeight','normal');
end

return
