function example_use
clear;clc;close all;

%initialize rates
beta = 1.1;
gamma = beta;
kinit= 2.1;
bs=10;

%initialize simulation parameters
nCells = 1e4;
rng(floor(log(4815162342)));
S = [bs 0; -1 1; 0 -1];
nT = 4;
kpar = [kinit,beta,gamma];
Tmax=5/min(kpar);
tvec = linspace(0,Tmax,nT);
t_matrix = repmat(tvec,nCells,1);

%run Gillespie simulation
X = gg_200110_gillespie_geom_1(kpar,t_matrix,S,nCells);


%initialize approximation parameters
N_approx_taylor = 8;
N_approx_laurent = 5;

dim = nT-1;
nrow=3;

tvec(end) = inf;

%visualize the Gillespie simulation results
for i = 1:(nT-1)
    %visualize mature marginal
    figure(1)
    subplot(1,dim,i)
    histogram(X(:,i+1,2),'BinMethod','integers','Normalization','pdf',...
        'FaceColor',0.5*[1 1 1],'EdgeColor','none'); hold on;
    title(sprintf('t = %.3f',tvec(i+1)),'FontWeight','Normal');
    xlabel('mRNA copy number');
    ylabel('Probability');
    
    %visualize nascent marginal
    figure(2)
    subplot(1,dim,i)
    histogram(X(:,i+1,1),'BinMethod','integers','Normalization','pdf',...
        'FaceColor',0.5*[1 1 1],'EdgeColor','none'); hold on;
    title(sprintf('t = %.3f',tvec(i+1)),'FontWeight','Normal');
    xlabel('pre-mRNA copy number');
    ylabel('Probability');
    
    %visualize joint distribution
    figure(3)
    subplot(nrow,dim,i)
    histogram2(X(:,i+1,1),X(:,i+1,2),'BinMethod','integers','Normalization','pdf',...
        'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none'); hold on;
    title(sprintf('t = %.3f',tvec(i+1)),'FontWeight','Normal');
    if i==1
        ylabel('Gillespie','FontWeight','bold');
    end
    if i==3
        yyaxis right
        set(gca,'ytick',[])
        ylabel('mRNA');
        set(gca,'ytick',[],'ycolor','k');
    end
    xlim([0,max(X(:,i+1,1))]);
    ylim([0,max(X(:,i+1,2))]); 
end

tvec(end)=inf;

%compute and visualize numerical integral results
for i = 1:(nT-1)
    M = max(X(:,i+1,1))+15;
    N = max(X(:,i+1,2))+15;
    
    %compute mature marginal
    Pa_marg =  gg_200228_numint_geom_tdep_3(kinit,bs,gamma,M,N,tvec(i+1),'mature');
    
    %visualize mature marginal
    figure(1)
    subplot(1,dim,i)
    plot(0:(N-1),Pa_marg,'b-','LineWidth',2);
    
    
    %compute nascent marginal
    Pa_marg =  gg_200228_numint_geom_tdep_3(kinit,bs,gamma,M,N,tvec(i+1),'nascent');
    
    %visualize nascent marginal
    figure(2)
    subplot(1,dim,i)
    plot(0:(M-1),Pa_marg,'b-','LineWidth',2);
    
    %compute joint distribution
    Pa =  gg_200228_numint_geom_tdep_3(kinit,bs,gamma,M,N,tvec(i+1),false);

    %visualize joint distribution
    figure(3)
    subplot(nrow,dim,i+dim);
    [x,y]=ndgrid((0:(M-1))-0.5,(0:(N-1))-0.5);
    h= pcolor(x,y,Pa);
    set(h, 'EdgeColor', 'none');
    xlim([0,max(X(:,i+1,1))]);
    ylim([0,max(X(:,i+1,2))]);
    if i==1
        ylabel('Numerical','FontWeight','bold');
    end
    if i==3
        yyaxis right
        set(gca,'ytick',[])
        ylabel('mRNA');
        set(gca,'ytick',[],'ycolor','k');
    end
end

%compute and visualize analytical integral results
for i = 1:(nT-1)
    M = max(X(:,i+1,1))+15;
    N = max(X(:,i+1,2))+15;
    
    %compute mature marginal
    Pa_marg =  gg_200325_analyt_geom_tdep_vec_31(kinit,bs,gamma,M,N,...
        tvec(i+1),'mature',N_approx_taylor,N_approx_laurent);
    
    %visualize mature marginal
    figure(1)
    subplot(1,dim,i)
    plot(0:(N-1),Pa_marg,'r--','LineWidth',2);
    ylim([min(ylim)/5,max(ylim)]);

    
    %compute nascent marginal
    Pa_marg =  gg_200325_analyt_geom_tdep_vec_31(kinit,bs,gamma,M,N,...
        tvec(i+1),'nascent',N_approx_taylor,N_approx_laurent);
    
    %visualize nascent marginal
    figure(2)
    subplot(1,dim,i)
    plot(0:(M-1),Pa_marg,'r--','LineWidth',2);
    ylim([min(ylim)/5,max(ylim)]);
        
    %compute joint distribution
    Pa =  gg_200325_analyt_geom_tdep_vec_31(kinit,bs,gamma,M,N,...
        tvec(i+1),false,N_approx_taylor,N_approx_laurent);
    
    figure(3)
    subplot(nrow,dim,i+2*dim);

    [x,y]=ndgrid((0:(M-1))-0.5,(0:(N-1))-0.5);
    h= pcolor(x,y,Pa);
    set(h, 'EdgeColor', 'none');
    xlim([0,max(X(:,i+1,1))]);
    ylim([0,max(X(:,i+1,2))]);
    
    %visualize joint distribution
    if i==1
        ylabel('Analytical','FontWeight','bold');
    end
    if i==3
        yyaxis right
        set(gca,'ytick',[])
        ylabel('mRNA');
        set(gca,'ytick',[],'ycolor','k');
    end
    xlabel('pre-mRNA');
end

set(figure(1),'Position',[254 561 824 319]);
set(figure(2),'Position',[254 561 824 319]);
set(figure(3),'Position',[402 316 654 606]);
return
