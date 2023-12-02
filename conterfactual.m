function res=conterfactual(y,A,B,MA,MB,H,MH,n,lags,M,N,T0,T1,a,nf)

TT0=(T0-63.25)*4;
TT1=(T1-63.25)*4;
a=((a-63.25))*4;
TT=TT1-TT0;
MM=M-N;
z=y;

YY=zeros(TT+3,n,MM);

BB1=zeros(3,3);
BB2=zeros(3,3);
CC=zeros(3,1);
AA=zeros(3,3);
HH=zeros(3,3);
se=zeros(3,1);

BBB1=zeros(3,3);
BBB2=zeros(3,3);
CCC=zeros(3,1);
ycf=zeros(TT+3,3);
ycf(1:2,:)=z(TT0-2:TT0-1,:);


for i=1:MM
    for t=TT0:TT1
        BB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);B(t,[16 17 18],i+N)];
        BB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);B(t,[19 20 21],i+N)];
        CC=(squeeze(B(t,[1 8 15],i+N)))';
        AA=tria(squeeze(A(t,:,i+N)));
        HH=diag(exp(squeeze(H(t,:,i+N))));
        se=HH\AA*(z(t,:)'-BB1*z(t-1,:)'-BB2*z(t-2,:)'-CC);
        
        if nf==5
            % heteroskedasticity
            %         BBB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);B(t,[16 17 18],i+N)];
            %         BBB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);B(t,[19 20 21],i+N)];
            %         CCC=(squeeze(B(t,[1 8 15],i+N)))';
            %         AAA=tria(squeeze(A(t,:,i+1+N)));
            %         HH=diag([exp(H(a(1),:,i+N+1))]);
            %         ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
            
            BBB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);B(t,[16 17 18],i+N)];
            BBB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);B(t,[19 20 21],i+N)];
            CCC=(squeeze(B(t,[1 8 15],i+N)))';
            AAA=tria(squeeze(A(t,:,i+1+N)));
            HH=diag([mean(exp(H(a(1)-7:a(1),:,i+N+1)),1)]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
        elseif nf==6
            
            BBB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);B(t,[16 17 18],i+N)];
            BBB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);B(t,[19 20 21],i+N)];
            CCC=(squeeze(B(t,[1 8 15],i+N)))';
            AAA=tria(squeeze(A(t,:,i+1+N)));
            HH=diag([exp(squeeze(H(t,1:2,i+N+1))) [mean(exp(H(a(1)-7:a(1),3,i+N+1)),1)]]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
        elseif nf==7
            
            BBB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);B(t,[16 17 18],i+N)];
            BBB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);B(t,[19 20 21],i+N)];
            CCC=(squeeze(B(t,[1 8 15],i+N)))';
            AAA=tria(squeeze(A(t,:,i+1+N)));
            HH=diag([[mean(exp(H(a(1)-7:a(1),1:2,i+N+1)),1)] exp(squeeze(H(t,3,i+N+1)))]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
            
        elseif nf==2
            % MP
            %         BBB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);B(a(1),[16 17 18],i+N)];
            %         BBB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);B(a(1),[19 20 21],i+N)];
            %         CCC=([squeeze(B(t,[1 8],i+N)) squeeze(B(a(1),[15],i+N))])';
            %         AAA=tria([squeeze(A(t,1,i+N+1)) squeeze(A(a(1),2:3,i+N+1))]);
            %         %HH=diag([exp(MH(a(1),:))]);
            %         ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
            BBB1=[B(t,[2 3 4],i+N);B(t,[9 10 11],i+N);mean(B(a(1)-7:a(1),[16 17 18],i+N),1)];
            BBB2=[B(t,[5 6 7],i+N);B(t,[12 13 14],i+N);mean(B(a(1)-7:a(1),[19 20 21],i+N),1)];
            CCC=([squeeze(B(t,[1 8],i+N)) squeeze(mean(B(a(1)-7:a(1),[15],i+N),1))])';
            AAA=tria([squeeze(A(t,1,i+N)) squeeze(mean(A(a(1)-7:a(1),2:3,i+N),1))]);
            %HH=diag([exp(MH(a(1),:))]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
            
        elseif nf==3
            % Private sector
            %         BBB1=[B(a(1),[2 3 4],i+N);B(a(1),[9 10 11],i+N);B(t,[16 17 18],i+N)];
            %         BBB2=[B(a(1),[5 6 7],i+N);B(a(1),[12 13 14],i+N);B(t,[19 20 21],i+N)];
            %         CCC=([squeeze(B(a(1),[1 8],i+N)) squeeze(B(t,[15],i+N))])';
            %         AAA=tria([squeeze(A(a(1),1,i+N+1)) squeeze(A(t,2:3,i+1+N))]);
            %         %HH=exp(squeeze(H(t,:,i+N+1)));
            %         %HH=diag([exp(MH(a(1),:))]);
            %         ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
            
            BBB1=[mean(B(a(1)-7:a(1),[2 3 4],i+N),1);mean(B(a(1)-7:a(1),[9 10 11],i+N),1);B(t,[16 17 18],i+N)];
            BBB2=[mean(B(a(1)-7:a(1),[5 6 7],i+N),1);mean(B(a(1)-7:a(1),[12 13 14],i+N),1);B(t,[19 20 21],i+N)];
            CCC=([squeeze(mean(B(a(1)-7:a(1),[1 8],i+N),1)) squeeze(B(t,[15],i+N))])';
            AAA=tria([squeeze(mean(A(a(1)-7:a(1),1,i+N+1),1)) squeeze(A(t,2:3,i+1+N))]);
            %HH=exp(squeeze(H(t,:,i+N+1)));
            %HH=diag([exp(MH(a(1),:))]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
        elseif nf==4
            % MP and private sector
            
            BBB1=[mean(B(a(1)-7:a(1),[2 3 4],i+N),1);mean(B(a(1)-7:a(1),[9 10 11],i+N),1);mean(B(a(1)-7:a(1),[16 17 18],i+N),1)];
            BBB2=[mean(B(a(1)-7:a(1),[5 6 7],i+N),1);mean(B(a(1)-7:a(1),[12 13 14],i+N),1);mean(B(a(1)-7:a(1),[19 20 21],i+N),1)];
            CCC=([squeeze(mean(B(a(1)-7:a(1),[1 8],i+N),1)) squeeze(mean(B(a(1)-7:a(1),[15],i+N),1))])';
            AAA=tria([squeeze(mean(A(a(1)-7:a(1),1,i+N+1),1)) squeeze(mean(A(a(1)-7:a(1),2:3,i+1+N),1))]);
            %HH=exp(squeeze(H(t,:,i+N+1)));
            %HH=diag([exp(MH(a(1),:))]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
        elseif nf==1
            % Everything
            %         BBB1=[B(a(1),[2 3 4],i+N);B(a(1),[9 10 11],i+N);B(a(1),[16 17 18],i+N)];
            %         BBB2=[B(a(1),[5 6 7],i+N);B(a(1),[12 13 14],i+N);B(a(1),[19 20 21],i+N)];
            %         CCC=([squeeze(B(a(1),[1 8],i+N)) B(a(1),15,i+N)])';
            %         AAA=tria([squeeze(A(a(1),1,i+N+1)) A(a(1),2:3,i+N+1)]);
            %         %HH=exp(squeeze(H(t,:,i+N+1)));
            %         HH=diag([exp(H(a(1),:,i+N+1))]);
            %         ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
            
            BBB1=[mean(B(a(1)-7:a(1),[2 3 4],i+N),1);mean(B(a(1)-7:a(1),[9 10 11],i+N),1);mean(B(a(1)-7:a(1),[16 17 18],i+N),1)];
            BBB2=[mean(B(a(1)-7:a(1),[5 6 7],i+N),1);mean(B(a(1)-7:a(1),[12 13 14],i+N),1);mean(B(a(1)-7:a(1),[19 20 21],i+N),1)];
            CCC=([squeeze(mean(B(a(1)-7:a(1),[1 8],i+N),1)) mean(B(a(1)-7:a(1),15,i+N),1)])';
            AAA=tria([squeeze(mean(A(a(1)-7:a(1),1,i+N+1),1)) mean(A(a(1)-7:a(1),2:3,i+N+1),1)]);
            %HH=exp(squeeze(H(t,:,i+N+1)));
            HH=diag([mean(exp(H(a(1)-7:a(1),:,i+N+1)),1)]);
            ycf(t-TT0+1+2,:)=(BBB1*ycf(t-1-TT0+1+2,:)'+BBB2*ycf(t-2-TT0+1+2,:)'+CCC+AAA\HH*se)';
        end
        
    end
    YY(:,:,i)=ycf(:,:);
end

res.YY=YY;
res.YYY=mean(YY,3);
res.YYYY=median(YY,3);
% plot([T0-.5:.25:T1],[res.YYY z(TT0-2:TT1,:)])
% std1=std(z(TT0-2:TT1,1));
% std2=std(z(TT0-2:TT1,2));
% stdc1=std(res.YY(:,1,:),1);
% stdc2=std(res.YY(:,2,:),1);
% portion1=1-stdc1./std1;
% portion2=1-stdc2./std2;
% res.medstd1=median(portion1);
% res.medstd2=median(portion2);
% res.meanstd1=mean(portion1);
% res.meanstd2=mean(portion2);
%
% mean1=mean(z(TT0:TT1,1));
% mean2=mean(z(TT0:TT1,2));
% meanc1=mean(res.YY(3:end,1,:),1);
% meanc2=mean(res.YY(3:end,2,:),1);
% hist(mean1-meanc1,20);
% figure
% hist(mean2-meanc2,20);
%
% res.m1=median(mean1-meanc1);
% res.sL1=prctile(mean1-meanc1,16);
% res.sU1=prctile(mean1-meanc1,84);
% res.me2=mean(mean2-meanc2);
% res.m2=median(mean2-meanc2);
% res.sL2=prctile(mean2-meanc2,16);
% res.sU2=prctile(mean2-meanc2,84);



% figure
% subplot(2,1,1); hist(portion1,20);
% subplot(2,1,2); hist(portion2,20);


CIL=prctile(squeeze(res.YY(:,1,:))',16);
CIU=prctile(squeeze(res.YY(:,1,:))',84);
figure
subplot(2,1,1);
plot([T0-.5:.25:T1],z(TT0-2:TT1,1),'LineWidth',2); grid; hold on;
plot([T0-.5:.25:T1],CIL','--r'); hold on;
plot([T0-.5:.25:T1],CIU','--r'); hold on;
plot([T0-.5:.25:T1],res.YYY(:,1),'k');
axis([70 88 0 12])
title('(a)')



CIL=prctile(squeeze(res.YY(:,2,:))',16);
CIU=prctile(squeeze(res.YY(:,2,:))',84);
subplot(2,1,2);
g1=plot([T0-.5:.25:T1],z(TT0-2:TT1,2),'LineWidth',2); grid; hold on;
g2=plot([T0-.5:.25:T1],CIL','--r'); hold on;
g3=plot([T0-.5:.25:T1],CIU','--r'); hold on;
g4=plot([T0-.5:.25:T1],res.YYY(:,2),'k');
axis([70 88 0 12])
title('(b)')
legend([g1 g2 g4],'actual','error bands','counterfactual')


% CIL=prctile(squeeze(res.YY(:,3,:))',16);
% CIU=prctile(squeeze(res.YY(:,3,:))',84);
% subplot(3,1,3);
% plot([T0-.5:.25:T1],z(TT0-2:TT1,3)); grid; hold on;
% plot([T0-.5:.25:T1],CIL',':r'); hold on;
% plot([T0-.5:.25:T1],CIU',':r'); hold on;
% plot([T0-.5:.25:T1],res.YYY(:,3),'c');




