function x=LASSO_fast(A,y,lambda)
[m,n]=size(A);
f=zeros(n,1);
x=zeros(n,1);
z=zeros(n,1);

mu=.5;

R=A'*A+eye(n)/mu;

iter_num=100;
sav=zeros(1,iter_num);
v=A'*y;

for iter=1:iter_num
    x=R\(v-f+z/mu);
    
    z=x+mu*f;
    temp=zeros(n,1);
    temp(z~=0)=z(z~=0)./abs(z(z~=0));
    r=abs(z)-lambda*mu;
    r(r<0)=0;
    z=r.*temp;
    
    f=f+mu*(x-z);
    sav(iter)=.5*norm(y-A*x)^2+lambda*norm(x,1);
end
