function impulse_resp_diff(AA,BB,HH,n,lags,M,N,dates,irf_dates)
% computes the difference between the impulse response in two different dates
% with 68% CI
% the dates are in the vector irf_dates and looked up using vector dates.
% these IR are computed NOT fixing the coefficients of the non-policy block to the 
% value of the first date (actually a mean over the previous year)


T=size(BB,1);
d=size(BB,2);
% Find time periods of selected irf_dates
a = arrayfun(@(d) find(dates==d), irf_dates);
[junk k]=size(a);
MM=M-N;

% Labels for graphs 
t_label = arrayfun(@(d) string(d, "yyyy:Q"), irf_dates); 
  
% Compute IRFs
H=diag(mean((HH(:,:)),1));
smat=zeros(3,3);
B=zeros(3,3,2);
response=zeros(k,n,n,20,MM);
for i=1:MM
for t=1:k
    smat=(tria([mean(squeeze(AA(a(t):a(t),1,i+N)),1);mean(squeeze(AA(a(t):a(t),2:3,i+N)),1)'])\H)';  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B(:,:,1)=[mean(BB(a(t):a(t),[2 3 4],i+N),1);mean(BB(a(t):a(t),[9 10 11],i+N),1);mean(BB(a(t):a(t),[16 17 18],i+N),1)];
    B(:,:,2)=[mean(BB(a(t):a(t),[5 6 7],i+N),1);mean(BB(a(t):a(t),[12 13 14],i+N),1);mean(BB(a(t):a(t),[19 20 21],i+N),1)];
    response(t,:,:,:,i)=impulsdtrf(B(:,:,:),smat,20);
   
end
end


for i=1:1 % Variable shocked (inflation) 
    for j=3:3 % Shock variable (interest rate, MP shock) 
        medp=prctile(squeeze(response(1,i,j,:,:))',50); 
        subplot(2,2,1); plot(medp,'-o'); grid; ; set(gca,'xtick',[0 5 10 15 20]); hold on; %axis([0 20 -.1 .1]);
        title('(a) Impulse responses of inflation')
        
        CILp=prctile(squeeze(response(1,i,j,:,:))'-squeeze(response(2,i,j,:,:))',16);
        CIUp=prctile(squeeze(response(1,i,j,:,:))'-squeeze(response(2,i,j,:,:))',84);
        medp=prctile(squeeze(response(1,i,j,:,:))'-squeeze(response(2,i,j,:,:))',50); 
        subplot(2,2,2); plot(CILp','--r'); grid; ; set(gca,'xtick',[0 5 10 15 20]); hold on; axis([0 20 -.1 .1]);
        plot(CIUp','--r'); plot(medp);
        title(sprintf('(b) Diff. %s and %s', t_label(1), t_label(2)))
        
        medp=prctile(squeeze(response(2,i,j,:,:))',50); 
        subplot(2,2,1); plot(medp,'-+'); grid; ; set(gca,'xtick',[0 5 10 15 20]); hold on; %axis([0 20 -.1 .1]); 
        CILp=prctile(squeeze(response(1,i,j,:,:))'-squeeze(response(3,i,j,:,:))',16);
        CIUp=prctile(squeeze(response(1,i,j,:,:))'-squeeze(response(3,i,j,:,:))',84);
        medp=prctile(squeeze(response(1,i,j,:,:))'-squeeze(response(3,i,j,:,:))',50); 
        subplot(2,2,3); plot(CILp','--r'); grid; ; set(gca,'xtick',[0 5 10 15 20]); hold on; axis([0 20 -.1 .1]); 
        plot(CIUp','--r'); plot(medp);
        title(sprintf('(c) Diff. %s and %s', t_label(1), t_label(3)))
        
        medp=prctile(squeeze(response(3,i,j,:,:))',50); 
        subplot(2,2,1); plot(medp,'-x'); grid; ; set(gca,'xtick',[0 5 10 15 20]); hold on; %axis([0 20 -.1 .1]); 
        legend(t_label{:})
        CILp=prctile(squeeze(response(2,i,j,:,:))'-squeeze(response(3,i,j,:,:))',16);
        CIUp=prctile(squeeze(response(2,i,j,:,:))'-squeeze(response(3,i,j,:,:))',84);
        medp=prctile(squeeze(response(2,i,j,:,:))'-squeeze(response(3,i,j,:,:))',50); 
        subplot(2,2,4); plot(CILp','--r'); grid; ; set(gca,'xtick',[0 5 10 15 20]); hold on; axis([0 20 -.1 .1]); 
        plot(CIUp','--r'); plot(medp);
        title(sprintf('(d) Diff. %s and %s', t_label(2), t_label(3)))
        
    end
end
% res=[];

