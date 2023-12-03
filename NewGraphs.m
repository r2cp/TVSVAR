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
time=[1963.5:.25:2001.5];


%% Figure 1: Posterior mean, 16-th and 84-th percentiles of the standard deviation of residuals of equations

figure(1)
set(gcf,"Position", [0 0 800 800]);
subplot(3,1,1);  
    plot(time,rsmall.HH(:,1), '-b'); hold on; title('(a) Standard deviation of residuals of the inflation equation')
    plot(time,prctile(exp(squeeze(rsmall.H(:,1,rsmall.N:end)))',16),'--r')
    plot(time,prctile(exp(squeeze(rsmall.H(:,1,rsmall.N:end)))',84),'--r')
    axis([1962 2002 0 .8])
    grid("on");
    
subplot(3,1,2); grid("on"); 
    plot(time,rsmall.HH(:,2), '-b'); hold on; title('(b) Standard deviation of residuals of the unemployment equation')
    plot(time,prctile(exp(squeeze(rsmall.H(:,2,rsmall.N:end)))',16),'--r')
    plot(time,prctile(exp(squeeze(rsmall.H(:,2,rsmall.N:end)))',84),'--r')
    axis([1962 2002 0 1.2])
    grid("on");
    
subplot(3,1,3); grid("on");
    plot(time,rsmall.HH(:,3), '-b'); hold on; title('(c) Standard deviation of residuals of interest rate equation (MP shocks)')
    plot(time,prctile(exp(squeeze(rsmall.H(:,3,rsmall.N:end)))',16),'--r')
    plot(time,prctile(exp(squeeze(rsmall.H(:,3,rsmall.N:end)))',84),'--r')
    axis([1962 2002 0 5])
    grid("on");


%% Figure 2: Impulse responses of inflation to monetary policy shocks 

% Periods to use
a=[75 81.5 96];
figure(2);
set(gcf,"Position", [0 0 800 800]);
impulse_resp_diff_infl(rsmall.A,rsmall.B,rsmall.HH,3,2,rsmall.M,rsmall.N, a);


%% Figure 3: Impulse responses of unemployment to monetary policy shocks 

figure(3)  
set(gcf,"Position", [0 0 800 800]);
impulse_resp_diff_unem(rsmall.A,rsmall.B,rsmall.HH,3,2,rsmall.M,rsmall.N,a);


%% Figure 4: Interest rate response to a 1% permanent increase of inflation 
% with 16th and 84th percentiles.

% computes the dinamic response of interest rate to a 1% shock to inflation
res = imp_strange(rsmall.A,rsmall.B,3,rsmall.M,rsmall.N,1);
axis([1960 2005 -.2 2.5])

%% Other figures... 
res=imp_strange(rsmall.A,rsmall.B,3,rsmall.M,rsmall.N,2);
axis([1960 2005 -2 0])

res=conterfactual(y(43:end,:),rsmall.A,rsmall.B,rsmall.AA,rsmall.BB,rsmall.H,rsmall.HH,3,2,rsmall.M,rsmall.N,70,87,92.75,2);


