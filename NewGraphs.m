% N=20001;
% r.BB=mean(r.B(:,:,N:end),3);
% r.AA=mean(r.A(:,:,N:end),3);
% r.HH=mean(exp(r.H(:,:,N:end)),3);
% time=[1963.5:.25:2001.5];


rsmall.B=r.B(:,:,[int:int:end]);
rsmall.A=r.A(:,:,[int:int:end]);
rsmall.H=r.H(:,:,[int:int:end]);

rsmall.M=size(rsmall.B,3);
rsmall.N=round(N/int);

rsmall.BB=mean(rsmall.B(:,:,rsmall.N:end),3);
rsmall.AA=mean(rsmall.A(:,:,rsmall.N:end),3);
rsmall.HH=mean(exp(rsmall.H(:,:,rsmall.N:end)),3);
time=[1963.5:.25:2001.5];





figure(1)
subplot(3,1,1); plot(time,rsmall.HH(:,1)); hold on; title('(a)')
                plot(time,prctile(exp(squeeze(rsmall.H(:,1,rsmall.N:end)))',16),'--r')
                plot(time,prctile(exp(squeeze(rsmall.H(:,1,rsmall.N:end)))',84),'--r')
                axis([1962 2002 0 .8])
                
subplot(3,1,2); plot(time,rsmall.HH(:,2)); hold on; title('(b)')
                plot(time,prctile(exp(squeeze(rsmall.H(:,2,rsmall.N:end)))',16),'--r')
                plot(time,prctile(exp(squeeze(rsmall.H(:,2,rsmall.N:end)))',84),'--r')
                axis([1962 2002 0 1.2])
                
subplot(3,1,3); plot(time,rsmall.HH(:,3)); hold on; title('(c)')
                plot(time,prctile(exp(squeeze(rsmall.H(:,3,rsmall.N:end)))',16),'--r')
                plot(time,prctile(exp(squeeze(rsmall.H(:,3,rsmall.N:end)))',84),'--r')
                axis([1962 2002 0 5])
                
a=[75 81.5 96];
figure(2)  
impulse_resp_diff_infl(rsmall.A,rsmall.B,rsmall.HH,3,2,rsmall.M,rsmall.N,a);
figure(3)  
impulse_resp_diff_unem(rsmall.A,rsmall.B,rsmall.HH,3,2,rsmall.M,rsmall.N,a);

res=imp_strange(rsmall.A,rsmall.B,3,rsmall.M,rsmall.N,1);
axis([1960 2005 -.2 2.5])
res=imp_strange(rsmall.A,rsmall.B,3,rsmall.M,rsmall.N,2);
axis([1960 2005 -2 0])

res=conterfactual(y(43:end,:),rsmall.A,rsmall.B,rsmall.AA,rsmall.BB,rsmall.H,rsmall.HH,3,2,rsmall.M,rsmall.N,70,87,92.75,2);

