function supp_fig_5
clear;clc;close all;
BATCH = 12;
SCAN=24;

p1=[227,172,82]/255;
p2=[252,222,164]/255;
p3=[90,180,172]/255;
n=8;
map = NaN(n*2,3);
for i = 1:3
    map(1:n,i) = linspace(p1(i),p2(i),n);
    map((n+1):end,i) = linspace(p2(i),p3(i),n);
end

figure(1);

NN=24;
T = zeros(NN,4,2500);
N_state = zeros(NN,1);

for DATANUM = 1:NN
    D = load(sprintf('data/gg_200618_sim_%i_%i.mat',BATCH,DATANUM));
    F = load(sprintf('landscape/gg_200625_land_scan_%i_%i_%i_%i.mat',SCAN,DATANUM,7,1));

    N_state(DATANUM) = (max(F.X(:,1))+1)*(max(F.X(:,2))+1);
    T(DATANUM,:,:) = [F.T_numint_chf_only,F.T_numint,...
        F.T_analytint_chf_only,F.T_analytint]';
end

t_names = {'Quad sans ifft2','Quad with ifft2','Specfun sans ifft2','Specfun with ifft2'};


for i = 1:4
    subplot(2,2,i)
    TT = squeeze(T(:,i,:))./N_state;
    TTlog = log10(TT);
    histogram(TTlog,'Normalization','Probability',...
        'FaceColor',0.5*[1 1 1],'EdgeColor','none');
    xlabel('log10 s per state');
    title(t_names{i},'FontWeight','Normal');
    xlim([-5.5,-4]);
    fprintf('%s: mean runtime per state %.3e s.\n',t_names{i},mean(TT(:)));
end

return
