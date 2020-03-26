function landscape_computer(N_T,N_L)
%load target Gillespie synthetic data
D = load('data/gill_syn_data_1.mat');
%only look at equilibrated timestep
X = squeeze(D.X(:,end,:));

%initialize grid of parameter values for evaluating the PMFs
kg_lim = [-1,1];
b_lim = [.1,2];
n_kg = 50; n_b = n_kg;
[KG_,B_] = ndgrid(logspace(kg_lim(1),kg_lim(2),n_kg),...
            logspace(b_lim(1),b_lim(2),n_b)); KG = KG_(:); B = B_(:);
N = length(KG);

%initialize empty variables to store information
T_numint = NaN(N,1); %time for numerical integration method
T_analytint = NaN(N,1); %time for analytical integration method
[div_ks,div_emd,div_lik,... %K-S, E-M, and K-L divergences btwn PMF outputs
data_ks_1,data_ks_2,... %K-S distance between outputs and target
data_emd_1,data_emd_2,...  %E-M distance between outputs and target
data_lik_1,data_lik_2] ... $K-L divergence between outputs and target
    = deal(NaN(N,1)); 

%this is fairly slow so best parallelized on a server
parfor i = 1:N
    %convert target copy-number data into joint histogram
    MAX = max(X,[],1)+1;
    h_data_pdf = histcounts2(X(:,1),X(:,2),...
        'BinMethod','integers','normalization','pdf',...
        'XBinLimits',[-0.5,MAX(1)-0.5],'YBinLimits',[-0.5,MAX(2)-0.5]);
    h_data_cdf = histcounts2(X(:,1),X(:,2),...
        'BinMethod','integers','normalization','cdf',...
        'XBinLimits',[-0.5,MAX(1)-0.5],'YBinLimits',[-0.5,MAX(2)-0.5]);


    %compute the numerical integral 
    tic 
    Pa =  gg_200128_numint_geom_tdep_2(KG(i),B(i),1,MAX(1),MAX(2),...
        inf,false);
    T_numint(i) = toc;
    %compute resulting output CDF
    Pa_n_CDF = cumsum(cumsum(Pa,1),2);

    %compute divergences from target data
    data_ks_1(i) = max(max(abs(Pa_n_CDF-h_data_cdf)));
    data_emd_1(i) = sum(sum(abs(Pa_n_CDF-h_data_cdf)));
    %round to avoid likelihood computation blowup
    EPS =1e-12;
    Pa(Pa<EPS)=EPS;
    data_lik_1(i) = sum(sum(h_data_pdf .* log(h_data_pdf./Pa),...
        'omitnan'),'omitnan');
    
    %store for later comparison with analytical integral
    Pa_num = Pa;

    %compute the analytical integral approximation using orders passed into
    %function
    tic 
    Pa =  gg_200130_analyt_geom_tdep_vec_26(KG(i),B(i),1,MAX(1),MAX(2),...
        inf,false,N_T,N_L);
    T_analytint(i) = toc;
    %compute resulting output CDF
    Pa_an_CDF = cumsum(cumsum(Pa,1),2);

    %compute divergences from target data
    %note here: since the calculated PMF is not guaranteed to be a
    %probability distribution, and the rounding occurs //after// the KS and
    %EMD calculations, the KS distance may be greater than 1. In practice,
    %the PMF would be rounded to zero or eps. This is somewhat of a
    %worst-case scenario.
    data_ks_2(i) = max(max(abs(Pa_an_CDF-h_data_cdf)));
    data_emd_2(i) = sum(sum(abs(Pa_an_CDF-h_data_cdf)));
    Pa(Pa<EPS)=EPS;
    %round to avoid likelihood computation blowup
    data_lik_2(i) = sum(sum(h_data_pdf .* log(h_data_pdf./Pa),...
        'omitnan'),'omitnan');

    %compare with numerical integral
    div_ks(i) = max(max(abs(Pa_n_CDF-Pa_an_CDF)));
    div_emd(i) = sum(sum(abs(Pa_n_CDF-Pa_an_CDF)));
    div_lik(i) = sum(sum(Pa_num .* log(Pa_num./Pa),'omitnan'),'omitnan');

end

%save all variables onto hard drive.
save(sprintf('landscape/gg_200131_land_scan_1_%i_%i.mat',N_T,N_L))
return
