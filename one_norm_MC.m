function t=one_norm_MC(p,beta,lambda, k)

epsilon=lambda*p/beta;

trial=1000000;
Gamma=randn(1,trial);
X=randn(1,trial);

S=p*Gamma+X;

Temp=S;
Temp(S>epsilon)=S(S>epsilon)-epsilon;
Temp(S<-epsilon)=S(S<-epsilon)+epsilon;
Temp((S<epsilon)&(S>-epsilon))=0;

t1=sum(abs(Temp-X))/trial;


S=p*Gamma;

Temp=S;
Temp(S>epsilon)=S(S>epsilon)-epsilon;
Temp(S<-epsilon)=S(S<-epsilon)+epsilon;
Temp((S<epsilon)&(S>-epsilon))=0;


t2=sum(abs(Temp))/trial;

t=k*t1+(1-k)*t2;