function [H snew]=step_sv(e,Vlast,in_shat,Hlast)
% This code makes a draw from the conditional posterior of the indicator
% variables and the volatility states 
%     e = observed innovations of the model with stochastic volatility
%     Vlast = last draw of the variance of the innovations to
%             the log standard deviations (must be a matrix)
%     in_shat = vector with guess 
%     Hlast = last draw of the volatility states
% Notice that here the log volatilities are assumed to follow random walk processes.
% This can be easily modified.

% setting parameters for the mixture of normals to approximate the logX^2
pr=[.0073 .10556 .00002 .04395 .34001 .24566 .25750]';
m=[-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819]'-1.2704;
v2=[5.79596 2.61369 5.1795 .16735 .64009 .34023 1.26261]';

% transforming the ss in a linear one
[T,n]=size(e);
e=log(e.^2+.001);   % I use the offset constant, as in KSC

% drawing the indicator variables s for the mixture of normals
s=zeros(T,n);
for i=1:n
    Q=repmat(pr',T,1).*normpdf(repmat(e(:,i),1,7),2*repmat(Hlast(:,i),1,7)+repmat(m',T,1),repmat(sqrt(v2)',T,1));  
    Q=Q./repmat(sum(Q,2),1,7);
    F=(Q+lag(Q',1)'+lag(Q',2)'+lag(Q',3)'+lag(Q',4)'+lag(Q',5)'+lag(Q',6)');
    states=unidrnd(100000,T,1)/100000;
    ss=sum((states*ones(1,7)>=F)')'+1;
    ss1=ss>7;
    s(:,i)=ss-ss1;
end

% forward recursion
SHAT=zeros(T,n);
SIG=zeros(T,n,n);
sig=eye(n);
shat=in_shat;
for t=1:T
    [shat,sig]=kfilter(e(t,:)-m(s(t,:))',2,1,shat,sig,diag(v2(s(t,:))),Vlast);
    SHAT(t,:)=shat';
    SIG(t,:,:)=sig;
end

% backward recursion
H=zeros(T,n);
H(T,:)=mvnrnd(shat,sig/2+sig'/2,1);
for t=T-1:-1:1
    [btTp,StTp] = kback(SHAT(t,:),squeeze(SIG(t,:,:)),H(t+1,:),1,Vlast);
    H(t,:) = mvnrnd(btTp,StTp/2+StTp'/2,1);   
end