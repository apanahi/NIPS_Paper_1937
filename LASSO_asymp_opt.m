function [s, beta, p]=LASSO_asymp_opt(lambda,delta,gamma,sigma2)

count=20;

beta=1;

s1=inner_opt(beta,lambda,delta,gamma,sigma2);
s2=inner_opt(beta/2,lambda,delta,gamma,sigma2);

flag_found=false;

if s1>s2
    betaup=beta;
    flag=true;
    c=1;
    while flag
        beta=beta/2;
        s1=s2;
        s2=inner_opt(beta/2,lambda,delta,gamma,sigma2);
        if (s2>s1)
            betadown=beta/2;
            flag=false;
        elseif (c==count)
            beta=beta/2;
            s=s2;
            flag=false;
            flag_found=true;
        end
        c=c+1;
    end
else
    betadown=beta/2;
    flag=true;
    c=1;
    while flag
        beta=2*beta;
        s2=s1;
        s1=inner_opt(beta,lambda,delta,gamma,sigma2);
        if (s1>s2)
            betaup=beta;
            flag=false;
        elseif (c==count)
            s=s1;
            flag=false;
            flag_found=true;
        end
        c=c+1;
    end
end

if ~flag_found
    for k=1:count
        a=(2*betaup+betadown)/3;
        b=(betaup+2*betadown)/3;
        sa=inner_opt(a,lambda,delta,gamma,sigma2);
        sb=inner_opt(b,lambda,delta,gamma,sigma2);
        if sa>sb
            betaup=a;
        else
            betadown=b;
        end
    end
    beta=(betaup+betadown)/2;
    [s,p]=inner_opt(beta,lambda,delta,gamma,sigma2);
end
s=-s;