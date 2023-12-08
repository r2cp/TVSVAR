function r=tvsvar(y,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVAR WITH TIME VARYING PARAMETERS AND VOLATILITIES
% These codes implement the estimation algorithm of Primiceri(2005) and 
% Del Negro and Primiceri (2015)
%
% The general form of the model is: y(t)=X(t)'*B(t)+D(t)*exp(H(t))*eps(t)
%                                   B(t)=B(t-1)+errB(t)
%                                   H(t)=H(t-1)+errH(t)
%                                   A(t)=D(t)^-1    (D and A lower triang. with ones on the main diag)
%                                   A(t)=A(t-1)+errA(t) (for the non zero/one elements of A)
%                                   V([errB(t) errA(t) errH(t)])=block diagonal
%                                   V(eps(t))=I
%
%
% y:    data matrix
%       (the constant is automatically added to the VAR)
%
% lags: number of lags.
%
% T0:   number of observations in the pre-sample used for the prior calibration
%
% T0B:  degrees of freedom for the IW prior distribution of var(errB)
%
% T0A:  vector of length n-1, containing the degrees of freedom for the IW
%       prior distributions of the blocks of var(errA)
%
% T0H:  degrees of freedom for the IW prior distribution of var(errH)
%
% kB:   constant scaling the prior of var(errB)
% kA:   constant scaling the prior of var(errA)
% kH:   constant scaling the prior of var(errH)
%
% M:    total number of draws in the Gibbs sampling algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin ~= 10); error('Wrong # of arguments'); end


% =========================================================================
% PRELIMINARY STEPS
% putting the VAR in the form Y=ZB+e and adding the constant to the regressors
x=[];
for i=1:lags
    x = [x lag(y,i)];
end
x=[ones(length(x),1) x];
x=x(lags+1:end,:);
y=y(lags+1:end,:);
[T,n]=size(y);
x0=[x(1:T0,:)];
x=[x(T0+1:end,:)];
y0=y(1:T0,:);
y=y(T0+1:end,:);
Y0=reshape(y0,n*T0,1);
Y=reshape(y,n*(T-T0),1);
Z0=kron(eye(n),x0);
Z=kron(eye(n),x);

% defining the relevant dimensions
T=T-T0;         % length of estimation sample
k = n + lags*n^2;   % # of VAR coefficients
na = n*(n-1)/2;   % # of free coefficients in the matrix A
ns = k+na+n;      % # of shocks with time-varying variance

% storage matrices
B=zeros(T,  k,   M);     
A=zeros(T,  na,  M+1);  
H=zeros(T,  n,   M+1);
V=zeros(ns, ns, M+1);


% =========================================================================
% SETTING UP THE PRIORS
% time varying coefficients in B
priorB=olsblock(y0,x0);
Bbar=reshape(priorB.bhatols,k,1)';
resB=priorB.resols;
VBbar=kron(resB'*resB/length(resB),eye(k/n)/priorB.XX);


% time varying coefficients in A
auxM=reshape([1:n^2],n,n)'; auxM=tril(auxM,-1); auxV=reshape(auxM,n^2,1); auxV(auxV==0)=[]; auxV=auxV';
ZZ0 = kron(eye(n),resB); ZZ0=ZZ0(:,auxV);
priorA = ols1(reshape(resB, n*T0, 1), ZZ0);
resA = reshape(priorA.resols, T0, n);
Abar=priorA.bhatols';
VAbarhat=(ZZ0'*ZZ0)\ZZ0'*(kron((resA'*resA/length(resA)),eye(T0)))*ZZ0/(ZZ0'*ZZ0);
VAbar=zeros(na);
count=1;
for j=1:n-1
    VAbar(count:count+j-1,count:count+j-1)=VAbarhat(count:count+j-1,count:count+j-1);
    count=count+j;
end

% time varying coefficients in H
Rbar = (diag(resA'*resA / length(resA)));
Hbar = log(Rbar)/2;

% auxiliary vector to be used below
auxMnew=reshape([1:(n-1)^2],n-1,n-1)'; auxMnew=tril(auxMnew); 
auxVnew=reshape(auxMnew',(n-1)^2,1); auxVnew(auxVnew==0)=[]; auxVnew=auxVnew';


% =========================================================================
% INITIALIZATION OF THE ALGORITHM

% Drawing the initial V from the prior
V0B=VBbar*T0B*kB^2;
V(1:k,1:k,1)=iwishrnd(V0B,T0B);
count=1;
for j=1:n-1
    V0A=T0A(j)*VAbar(count:count+j-1,count:count+j-1)*kA^2;
    V(k+count:k+count+j-1,k+count:k+count+j-1,1)=iwishrnd(V0A,T0A(j));
    count=count+j;
end

% % For constant volatilities
% % Drawing the initial W = Var(errH)
% V0H = eye(n)*T0H*kH^2;
% V(k+na+1:end, k+na+1:end, 1) = iwishrnd(V0H,T0H);

% Volatilities constant through time 
V(k+na+1:end, k+na+1:end, 1) = zeros(n, n);

% Drawing the initial H
H(1,:,1) = mvnrnd(Hbar, eye(n), 1);
for t=2:T
    % Normal distribution with mean H(t-1,:,1) and variance V(k+na+1:end,k+na+1:end,1)
    % H(t,:,1) = mvnrnd(H(t-1,:,1), V(k+na+1:end,k+na+1:end,1), 1);

    % Volatilities constant through time
    H(t, :, 1) = H(t-1, :, 1); 
end

% drawing the initial A from the prior
A(1,:,1)=mvnrnd(Abar,4*(VAbar/2+VAbar'/2),1);
for t=2:T
    A(t,:,1)=mvnrnd(A(t-1,:,1),V(k+1:k+na,k+1:k+na,1),1);
end

% =========================================================================
% =========================================================================
% GIBBS SAMPLING ALGORITHM

for i=1:M
    if i==20*floor(.05*i);
    % if mod(i, floor(0.05*i)) == 0
        disp(i);
    end
    
    % =====================================================================
    % STEP 1: DRAWS OF B

    % Kalman filter
    SHAT=zeros(T,k);
    SIG=zeros(T,k,k);
    sig=4*VBbar;
    shat=Bbar';
    for t=1:T
        [shat,sig] = kfilter (y(t,:)',Z([0:n-1]*T+t,:),eye(k),shat,sig,...
            tria(A(t,:,i))\diag(exp(2*squeeze(H(t,:,i))))/(tria(A(t,:,i))'),squeeze(V(1:k,1:k,i)));
        SHAT(t,:)=shat';
        SIG(t,:,:)=sig;
    end
    
    % simulation smoother
    B(T,:,i)=mvnrnd(shat,sig/2+sig'/2,1);
    for t=T-1:-1:1
        [btTp,StTp]=kback(SHAT(t,:),squeeze(SIG(t,:,:)),squeeze(B(t+1,:,i)),eye(k),squeeze(V(1:k,1:k,i)));
        B(t,:,i)=mvnrnd(btTp,StTp/2+StTp'/2,1);
    end


    % =====================================================================
    % STEP 2: DRAWS OF A

    Yhat = (Y-sum((Z.*repmat(squeeze(B(:,:,i)),n,1))')');
    yhat = reshape(Yhat,T,n);
    
    % Kalman filter
    SHAT=zeros(T,na);
    SIG=zeros(T,na,na);
    sig=4*VAbar;
    shat=Abar';
    for t=1:T
        yreg=kron(eye(n-1),yhat(t,1:end-1));
        yreg=yreg(:,auxVnew);

        [shat,sig] = kfilter (yhat(t,2:end)',yreg,eye(na),shat,sig,diag(exp(2*squeeze(H(t,2:end,i)))),squeeze(V(k+1:k+na,k+1:k+na,i)));
        SHAT(t,:)=shat';
        SIG(t,:,:)=sig;
    end

    % simulation smoother
    A(T,:,i+1) = mvnrnd(shat,sig/2+sig'/2,1);
    for t=T-1:-1:1
        [btTp,StTp] = kback(SHAT(t,:),squeeze(SIG(t,:,:)),squeeze(A(t+1,:,i+1)),eye(na),squeeze(V(k+1:k+na,k+1:k+na,i)));
        A(t,:,i+1) = mvnrnd(btTp,StTp/2+StTp'/2,1);
    end

    
    % =====================================================================
    % STEP 3: DRAWS OF the indicator variables s and H
    ystar = zeros(T,n);
    for t=1:T
        ystar(t,:) = (tria(A(t,:,i+1)) * yhat(t,:)')';
    end
  
    % Sigma matrix is obtained for all periods
    % [H(:,:,i+1)] = step_sv(ystar, squeeze(V(k+na+1:end,k+na+1:end,i)), Hbar, squeeze(H(:,:,i)));
    % [H(:,:,i+1)] = step_sv_mcmc(ystar,squeeze(V(k+na+1:end,k+na+1:end,i)), Hbar, squeeze(H(:,:,i))); % use this line instead of the previous one to implement the exact algorithm (algorithm 3)    
    
    % Diagonal entries of Sigma matrix are constant under W = var(errH) = 0 (nxn)
    % W = V(k+na+1:end, k+na+1:end, i)
    Haux = step_sv(ystar, squeeze(V(k+na+1:end, k+na+1:end, i)), Hbar, squeeze(H(:,:,i)));
    H(:,:,i+1) = Haux;
    
    % =====================================================================
    % STEP 4: DRAWS OF V
    
    % Drawing var(errB)
    errB=[B(2:end,:,i)-B(1:end-1,:,i)];
    V1B=errB'*errB+V0B;
    V(1:k,1:k,i+1)=iwishrnd(V1B,T-1+T0B);

    % Drawing var(errA)
    count=1;
    for j=1:n-1
        errA = A(2:end, count:count+j-1, i+1) - A(1:end-1, count:count+j-1, i+1);
        V1A = errA' * errA + T0A(j) * VAbar(count:count+j-1,count:count+j-1) * kA^2;
        % Fill S = Var(A)
        V(k+count:k+count+j-1, k+count:k+count+j-1, i+1) = iwishrnd(V1A, T-1+T0A(j));
        count = count+j;
    end
    
    % Drawing var(errH)
    % errH = H(2:end,:,i+1) - H(1:end-1,:,i+1);
    % V1H = errH' * errH + (kH^2) * eye(n) * T0H;
    % V(k+na+1:end, k+na+1:end, i+1) = iwishrnd(V1H, T-1+T0H);

    % Removing the process of noise W from V 
    V(k+na+1:end, k+na+1:end, i+1) = zeros(n, n); 

end

r.B=B;
r.V=V;
r.A=A;
r.H=H;