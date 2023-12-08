close all;
clear all;
clc; 

load data          % Original paper data through 2001
load recentdata    % Data through 2019:4

%% TVSVAR Model 

y=recentdata;
lags=2;
T0=40;      % number of observations in the pre-sample used for the prior calibration
T0B=40;     % degrees of freedom for the IW prior distribution of var(errB)
T0A=[2 3];  % vector of length n-1, containing the degrees of freedom for the IW prior distributions of the blocks of S =   var(errA)
T0H=4;      % degrees of freedom for the IW prior distribution of var(errH)
kB=.01;     % constant scaling the prior of Q = var(errB)
kA=.1;      % constant scaling the prior of S = var(errA)
kH=.01;     % constant scaling the prior of W = var(errH)
M=10000;      % total number of draws in the Gibbs sampling algorithm

% Run the Time-Varying SVAR MCMC estimation
r = tvsvar(y,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);

%% Load data for graphs

% Replication of part (a), with data through 2001
% load workspace_data_10000.mat

% Part (b), with data through 2019Q4
% load workspace_recentdata_10000.mat

% Part (c), with data through 2019Q4 - Constant volatilities
% load workspace_recentdata_constvol_10000.mat

%% Generate graphs
close all;

int=1;      % thinning parameter (choose 1 if want to use all draws for constructing 
            % the final graphs; choose J if want to use one every J draws for 
            % constructing the final graphs) 
N=2000;     % # of draws to burn

t1 = datetime(1953,1,1)+ (T0 + lags)*calmonths(3); 
t2 = datetime(2019,12,1);
dates = t1:calmonths(3):t2;

NewGraphs;