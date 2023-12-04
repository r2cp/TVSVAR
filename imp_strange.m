function res=imp_strange(A,B,n,M,N,k)
% computes the dynamic response of interest rate to a 1% shock of variable
% k
% note: k denotes the variable of the permanent shock: 
% k=1 -> inflation, k=2 -> unemployment

a=[1 10 20 60];
MM=M-N;
T=size(B,1);
B1=zeros(T,n,n);
B2=zeros(T,n,n);
I=zeros(T,60);
II=zeros(T,4,MM);
for i=1:MM
    AA=squeeze(A(:,:,N+i));
    BB=squeeze(B(:,:,N+i));
    for t=1:T
        B1(t,:,:)=tria(AA(t,:))*[BB(t,[2 3 4]);BB(t,[9 10 11]);BB(t,[16 17 18])];
        B2(t,:,:)=tria(AA(t,:))*[BB(t,[5 6 7]);BB(t,[12 13 14]);BB(t,[19 20 21])];
    end
    
    
    I(:,1)=AA(:,k+1);
    I(:,2)=squeeze(B1(:,3,3)).*I(:,1)+squeeze(B1(:,3,k))+AA(:,k+1);
    I(:,3)=squeeze(B1(:,3,3)).*I(:,2)+squeeze(B2(:,3,3)).*I(:,1)+squeeze(B1(:,3,k))+AA(:,k+1)+squeeze(B2(:,3,k));
    for j=4:60
        I(:,j)=squeeze(B1(:,3,3)).*I(:,j-1)+squeeze(B2(:,3,3)).*I(:,j-2)+squeeze(B1(:,3,k))+AA(:,k+1)+squeeze(B2(:,3,k));
    end
    II(:,:,i)=I(:,a);
end
res.II=II;

figure
set(gcf,"Position", [0 0 800 800]);
dates = 1963.5:.25:2019.75; % 1963Q3 - 2019Q4

for j=1:4
    
    CILp=prctile(squeeze(II(:,j,:))',16);
    CIUp=prctile(squeeze(II(:,j,:))',84);
    medp=prctile(squeeze(II(:,j,:))',50);
    %figure
    %         if j==4
    %         figure(5);
    %         plot([1963.25:.25:2001.75]',[medp']); grid;
    %         hold on
    %         plot([1963.25:.25:2001.75]',CILp','--r');
    %         plot([1963.25:.25:2001.75]',CIUp','--r');
    %         if k==1
    %             axis([1962 2002 -1 3.5])
    %         else
    %             axis([1962 2002 -3 1])
    %         end
    %     else
    subplot(2,2,j);
    plot(dates', [medp']); grid;
    hold on
    plot(dates',CILp','--r');
    plot(dates',CIUp','--r');
    if k==1 % For permanent shock of inflation 
        if j ~= 4
            % lims for panels (a)-(c)
            axis([1962 2020.25 -1 3.5])
        else
            % lims for panel (d) plot
            axis([1962 2020.25 -1 4.2])
        end
    else % For permanent shock of unemployment
        axis([1962 2020.25 -3 1])
        if j == 4 
            axis([1962 2020.25 -3.5 1])
        end
    end
    %end
end
subplot(2,2,1); title('(a) Simultaneous response');
subplot(2,2,2); title('(b) Response after 10 quarters');
subplot(2,2,3); title('(c) Response after 20 quarters');
subplot(2,2,4); title('(d) Response after 60 quarters');

% Figure 5: with just the median responses
figure
set(gcf,"Position", [0 0 800 800]);
me=median(II(:,:,:),3);
plot(dates',me(:,1),'--'); grid; hold on
plot(dates',me(:,2),':');  hold on
plot(dates',me(:,3),'-.');  hold on
plot(dates',me(:,4));
legend('Response after 0 quarters','Response after 10 quarters','Response after 20 quarters','Response after 60 quarters')
if k == 1 % lims for inflation
    axis([1960 2020 -.2 3.5])
else % lims for unemployment
    axis([1960 2020 -2.2 0])
end