%1/6/2020
%unif
function X_mesh = ...
    gg_200106_gillespie_unif_1(k,t_matrix,S,nCells)

num_t_pts = size(t_matrix,2);

X_mesh = NaN(nCells,num_t_pts,2); 

t = zeros(nCells,1); 
tindex = ones(nCells,1);

%initialize state: integer unspliced, integer spliced 
X = zeros(nCells,2);

%initialize list of cells that are being simulated
simindices = 1:nCells;
activecells = true(nCells,1);
%define offset time to randomize starting position of deterministic
%initialization rate
% offset = rand(nCells,1) * (1/kon + 1/koff);
% offset = zeros(nCells,1);

while any(activecells)
    mu = zeros(nCells,1);
    
    [t_upd,mu_upd] = rxn_calculator(...
        X(activecells,:),...
        t(activecells),...
        k,...
        sum(activecells));

    t(activecells) = t_upd;
    mu(activecells) = mu_upd;
    
    linindupdate = sub2ind(size(t_matrix),(1:length(tindex(activecells)))',...
        tindex(activecells));
    tvec_time = t_matrix(linindupdate);
    update = false(nCells,1);
    update(activecells) = t(activecells)>tvec_time;
%     end
    
    while any(update)
        tobeupdated = find(update);
        for i = 1:length(tobeupdated)
            X_mesh(simindices(tobeupdated(i)),tindex(tobeupdated(i)),:) = ...
                double(X(tobeupdated(i),:));
        end
        tindex = tindex+update;
        ended_in_update = tindex(update)>num_t_pts;

        if any(ended_in_update)
            ended = tobeupdated(ended_in_update);
            
            activecells(ended) = false;
            mu(ended) = 0;

            if ~any(activecells),break;end
        end
        
        linindupdate = sub2ind(size(t_matrix),(1:length(tindex(activecells)))',...
            tindex(activecells));
        tvec_time = t_matrix(linindupdate);
        update = false(nCells,1);
        update(activecells) = t(activecells)>tvec_time;

    end
    
    z_ = find(activecells);
    not_burst = mu(z_) > 1;
    burst = mu(z_) == 1;
%     X(activecells,:) = X(activecells,:) + S(mu(activecells),:);
    X(z_(not_burst),:) = X(z_(not_burst),:) + S(mu(z_(not_burst)),:);
    bs = geornd(1/(1+S(1,1)),[sum(burst),1]);
    
%     [bs zeros(sum(burst),1)]
    
    X(z_(burst),:) = X(z_(burst),:) + [bs zeros(sum(burst),1)];
    
    
end

return


function [t,mu] = rxn_calculator(X,t,k,nCells)
nRxn = 3;

a = zeros(nCells,nRxn);

% a is propensity matrix
% reactions:
% production
% conversion
% death

kinit = k(1);
beta = k(2);
gamma = k(3);

a(:,1) = kinit;
a(:,2) = beta * X(:,1);
a(:,3) = gamma * X(:,2);

a0 = sum(a,2);
t = t + log(1./rand(nCells,1))./a0;
r2ao = a0.*rand(nCells,1);
mu = sum(repmat(r2ao,1,nRxn+1) >= cumsum([zeros(nCells,1),a],2),2);
return
