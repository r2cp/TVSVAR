function [btT,StT]=kback(btt,Stt,bnT,A,omega)
%[btT StT]=ksmooth(btt,Stt,bnT,SnT,A,omega)
% Smoothing recursion.  State evolution equation is
%    bn=A*bt+e,  Var(e)=omega
%    bt|t ~ N(btt,Stt) -- from Kalman Filter
%    bn|T ~ N(bnT,SnT) -- distribution of bn given full sample. From 
%                         KF if n=T, otherwise from this recursion
%    bt|T ~ N(btT,StT)
AS=A*Stt;
G=AS*A'+omega;
SAGI=AS'/G;
btT=(SAGI*(bnT'-A*btt'))'+btt;
StT=Stt-SAGI*AS;%SAGI*SnT*SAGI';
