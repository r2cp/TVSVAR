close all;
clear all;

load data          % Original paper data through 2001
load recentdata    % Data through 2019:4

%% TVSVAR Model 
y=data;
lags=2;
T0=40;      % number of observations in the pre-sample used for the prior calibration
T0B=40;     % degrees of freedom for the IW prior distribution of var(errB)
T0A=[2 3];  % vector of length n-1, containing the degrees of freedom for the IW prior distributions of the blocks of var(errA)
T0H=4;      % degrees of freedom for the IW prior distribution of var(errH)
kB=.01;     % constant scaling the prior of Q = var(errB)
kA=.1;      % constant scaling the prior of S = var(errA)
kH=.01;     % constant scaling the prior of W = var(errH)
M=10000;      % total number of draws in the Gibbs sampling algorithm

% Run the Time-Varying SVAR MCMC estimation
r = tvsvar(y,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);

%% Graphs
close all;
% load workspace_data_10000.mat

int=1;     % thinning parameter (choose 1 if want to use all draws for constructing 
            % the final graphs; choose J if want to use one every J draws for 
            % constructing the final graphs) 
N=2000;
NewGraphs;