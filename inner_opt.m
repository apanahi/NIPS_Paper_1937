function [s, p]=inner_opt(beta,lambda,delta,gamma,sigma2)
count=20;

p=1;

s1=ExpProx(beta,p,lambda,delta,gamma,sigma2);
s2=ExpProx(beta,p/2,lambda,delta,gamma,sigma2);

flag_found=false;

if s1>s2
    pup=p;
    flag=true;
    c=1;
    while flag
        p=p/2;
        s1=s2;
        s2=ExpProx(beta,p/2,lambda,delta,gamma,sigma2);
        if (s2>s1)
            pdown=p/2;
            flag=false;
        elseif (c==count)
            p=p/2;
            s=s2;
            flag=false;
            flag_found=true;
        end
        c=c+1;
    end
else
    pdown=p/2;
    flag=true;
    c=1;
    while flag
        p=2*p;
        s2=s1;
        s1=ExpProx(beta,p,lambda,delta,gamma,sigma2);
        if (s1>s2)
            pup=p;
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
        a=(2*pup+pdown)/3;
        b=(pup+2*pdown)/3;
        sa=ExpProx(beta,a,lambda,delta,gamma,sigma2);
        sb=ExpProx(beta,b,lambda,delta,gamma,sigma2);
        if sa>sb
            pup=a;
        else
            pdown=b;
        end
    end
    p=(pup+pdown)/2;
    s=ExpProx(beta,p,lambda,delta,gamma,sigma2);    
end

s=-s;