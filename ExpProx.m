function s=ExpProx(beta,p,lambda,delta,gamma,sigma2)

q=beta/p;

sigma=p;

e=exp(-(lambda/sigma/q)^2/2);
Q=qfunc(lambda/sigma/q);
s1=2*(lambda*sigma*e/2/sqrt(2*pi)-Q*(lambda^2/2/q+q*sigma^2/2)+q*sigma^2/4);

sigma=sqrt(1+p^2);

e=exp(-(lambda/sigma/q)^2/2);
Q=qfunc(lambda/sigma/q);
s2=2*(lambda*sigma*e/2/sqrt(2*pi)-Q*(lambda^2/2/q+q*sigma^2/2)+q*sigma^2/4);

s=(1-delta)*s1+delta*s2+gamma*sigma2*beta/2/p+(gamma-1)*beta*p/2-gamma*beta^2/2;
