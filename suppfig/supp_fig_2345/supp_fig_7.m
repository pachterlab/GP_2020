function supp_fig_7
clear;clc;close all;
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
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
divs = zeros(NN,4,2500);
N_state = zeros(NN,1);

for DATANUM = 1:NN
    D = load(sprintf('data/gg_200618_sim_%i_%i.mat',BATCH,DATANUM));
    F = load(sprintf('landscape/gg_200625_land_scan_%i_%i_%i_%i.mat',SCAN,DATANUM,7,1));

    N_state(DATANUM) = (max(F.X(:,1))+1)*(max(F.X(:,2))+1);
    divs(DATANUM,:,:) = [F.div_lik,F.div_ks,...
        F.div_emd,F.div_chf]';
end

t_names = {'KL divergence','KS distance','log10 Wasserstein distance','log10 chf distance'};

N_state = repmat(N_state,[1,2500]);

for i = 1:4
    subplot(2,2,i)
    TT = squeeze(divs(:,i,:));
    if i>2
        TT = log10(TT);
    end
    scatter(log10(N_state(:)),TT(:),10,'k','filled','MarkerFaceAlpha',0.05)
    ylabel(t_names{i},'FontWeight','Normal');
    xlabel('log10 state space size');
end

return
