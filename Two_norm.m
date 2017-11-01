function t=Two_norm(p,beta,lambda, k)

epsilon=lambda*p/beta;

t1=2*qfunc(epsilon/p)*(p^2+epsilon^2)-2*exp(-epsilon^2/2/p^2)*epsilon*p/sqrt(2*pi);

t2=2*qfunc(epsilon/sqrt(1+p^2))*(p^2+epsilon^2-1)-2*exp(-epsilon^2/2/(1+p^2))*epsilon*sqrt(1+p^2)/sqrt(2*pi)+1;

t=t1*(1-k)+t2*k;