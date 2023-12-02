function M=tria(x);
m=length(x);
n=(1+sqrt(1+8*m))/2;
M=eye(n);
M(logical(triu(ones(n),1)))=-x;
M=M';