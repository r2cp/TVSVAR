% N=20001;
% r.BB=mean(r.B(:,:,N:end),3);
% r.AA=mean(r.A(:,:,N:end),3);
% r.HH=mean(exp(r.H(:,:,N:end)),3);
% time=[1963.5:.25:2001.5];


%% Obtaining a smaller subsample of the posterior for the graphs

% Subsample of the posterior parameters
rsmall.B=r.B(:,:,[int:int:end]);
rsmall.A=r.A(:,:,[int:int:end]);
rsmall.H=r.H(:,:,[int:int:end]);

rsmall.M=size(rsmall.B,3);
rsmall.N=round(N/int);

% Mean of the parameters and standard deviation (\sigma_{it})
rsmall.BB=mean(rsmall.B(:,:,rsmall.N:end),3);
rsmall.AA=mean(rsmall.A(:,:,rsmall.N:end),3);
rsmall.HH=mean(exp(rsmall.H(:,:,rsmall.N:end)),3);

% Time periods
% e.g. 1963.5 = 1963Q2, 2019.75 = 2019Q4
time=[1963.5:.25:2019.75];


%% Figure 1: Posterior mean, 16-th and 84-th percentiles of the standard deviation of residuals of equations

figure(1)
set(gcf,"Position", [0 0 800 800]);
subplot(3,1,1);  
    plot(time,rsmall.HH(:,1), '-b'); hold on; title('(a) Standard deviation of residuals of the inflation equation')
    plot(time,prctile(exp(squeeze(rsmall.H(:,1,rsmall.N:end)))',16),'--r')
    plot(time,prctile(exp(squeeze(rsmall.H(:,1,rsmall.N:end)))',84),'--r')
    axis([1962 2020 0 .8])
    grid("on");
    
subplot(3,1,2); grid("on"); 
    plot(time,rsmall.HH(:,2), '-b'); hold on; title('(b) Standard deviation of residuals of the unemployment equation')
    plot(time,prctile(exp(squeeze(rsmall.H(:,2,rsmall.N:end)))',16),'--r')
    plot(time,prctile(exp(squeeze(rsmall.H(:,2,rsmall.N:end)))',84),'--r')
    axis([1962 2020 0 1.2])
    grid("on");
    
subplot(3,1,3); grid("on");
    plot(time,rsmall.HH(:,3), '-b'); hold on; title('(c) Standard deviation of residuals of interest rate equation (MP shocks)')
    plot(time,prctile(exp(squeeze(rsmall.H(:,3,rsmall.N:end)))',16),'--r')
    plot(time,prctile(exp(squeeze(rsmall.H(:,3,rsmall.N:end)))',84),'--r')
    axis([1962 2020 0 5])
    grid("on");


%% Figure 2: Impulse responses of inflation to monetary policy shocks 

% Periods to use to show differences in IRFs
irf_dates = [datetime(1996,1,1), datetime(2010,1,1), datetime(2017,1,1)]; 

figure(2);
set(gcf,"Position", [0 0 800 800]);
impulse_resp_diff_infl(rsmall.A,rsmall.B,rsmall.HH,3,2,rsmall.M,rsmall.N, dates, irf_dates);


%% Figure 3: Impulse responses of unemployment to monetary policy shocks 

figure(3)  
set(gcf,"Position", [0 0 800 800]);
impulse_resp_diff_unem(rsmall.A,rsmall.B,rsmall.HH,3,2,rsmall.M,rsmall.N, dates, irf_dates);


%% Figure 4 and 5: Interest rate response to a 1% permanent increase of inflation 
% with 16th and 84th percentiles.

% computes the dinamic response of interest rate to a 1% shock of inflation
res = imp_strange(rsmall.A,rsmall.B,3,rsmall.M,rsmall.N, 1);

%{
%% Figures 6 and 7: Interest rate response to a 1% permanent increase of unemployment

% computes the dinamic response of interest rate to a 1% shock of
% unemployment
res = imp_strange(rsmall.A,rsmall.B,3,rsmall.M,rsmall.N, 2);

%% Figure 8: Counterfactual historical simulation 
% drawing the parameters of the monetary policy rule from their 1991-1992
% posterior and using it in the 1970s 
% 
% "Drawing the parameters of the policy rule in the 1970's from their
% posterior in 1991-1992, in order to see whether this would have made any
% difference."

% Still need to check and modify here for the extended sample...
res=conterfactual(y((T0 + lags + 1):end, :),rsmall.A,rsmall.B,rsmall.AA,rsmall.BB,rsmall.H,rsmall.HH, ... 
    3, 2, ... % Variables and lags
    ... rsmall.M,rsmall.N, ... % Draws and burn
    10000, 2000, ... 
    dates(27:end), ...  % dates
    70, ...             % T0 = 1970Q1
    119.75, ...         % T1 = 2019Q4
    ... 92.75, ...          % a (1992Q4, specific date to take parameters) 
    117.75, ...          % a (2017Q4, specific date to take parameters) 
    2 ...               % nf = scenario number?
); 


%}
