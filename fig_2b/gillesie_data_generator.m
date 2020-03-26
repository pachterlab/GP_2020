function gillesie_data_generator
clear;clc;close all;
%%%%%%%
% Set physiology parameters
beta = 1;
gamma = beta;
kinit= 2.5;
bs=19;

%set Gillespie simulation parameters
nCells = 1e4;
rng(floor(log(4815162342)));
S = [bs 0; -1 1; 0 -1];
nT = 4;
kpar = [kinit,beta,gamma];
Tmax=5/min(kpar);
tvec = linspace(0,Tmax,nT);
t_matrix = repmat(tvec,nCells,1);

%%%%%
% Run simulation
X = gg_200110_gillespie_geom_1(kpar,t_matrix,S,nCells);

%save output
save('data/gill_syn_data_1.mat','kpar','S','X')
return
