
delta=.2;
gamma=.5;
sigma2=1e-1;

n=200;
m=fix(gamma*n);

trial_num=1000;


c=0;

lambda_vec=.01:1:5.01;
savs=zeros(1,length(lambda_vec));
savM=zeros(1,length(lambda_vec));
savtwo=zeros(1,length(lambda_vec));
savone=zeros(1,length(lambda_vec));
savbeta=zeros(1,length(lambda_vec));
savp=zeros(1,length(lambda_vec));
valsav=zeros(1,length(lambda_vec));
valerror=zeros(1,length(lambda_vec));
varerror=zeros(1,length(lambda_vec));

%%% For nonsymetric bernouli %%%
pnb=.1;
anb=sqrt((1-pnb)/pnb);
bnb=-sqrt(pnb/(1-pnb));
%%%%%%

A=zeros(m,n);
for lambda=lambda_vec
    c=c+1
    [s, beta, p]=LASSO_asymp_opt(lambda,delta,gamma,sigma2);
    
    savs(c)=s;
    savbeta(c)=beta;
    savp(c)=p;
    savtwo(c)=Two_norm(p,beta,lambda, delta);
    savM(c)=(1-delta)*qfunc(lambda/beta)+ delta*qfunc(lambda*p/beta/sqrt(1+p^2));
    savone(c)=one_norm_MC(p,beta,lambda, delta);

    
    temp_val=0;
    temp_error=0;
    temp_error_sqr=0;
    for trial=1:trial_num
        Anb=(rand(m,n)<pnb);
        A(Anb==1)=anb/sqrt(m);
        A(Anb==0)=bnb/sqrt(m);
        
        %A=sqrt(3-2)*trnd(3,m,n)/sqrt(m)/sqrt(3);
        
        x0=randn(n,1);
        x0(rand(n,1)>delta)=0;
        
        y=A*x0+sqrt(sigma2)*randn(m,1);
        
        x=LASSO_fast(A,y,lambda);
        temp_val=temp_val+(.5*norm(y-A*x)^2+lambda*norm(x,1))/n;
        temp_error=temp_error+norm(x-x0)^2/n;
        temp_error_sqr=temp_error_sqr+(norm(x-x0)^2/n)^2;
    end
    valsav(c)=temp_val/trial_num;
    valerror(c)=temp_error/trial_num;
    varerror(c)=temp_error_sqr/trial_num-valerror(c)^2;
end

figure(1)
hold on
plot(lambda_vec,savtwo,lambda_vec,valerror,'--o',...
    'LineWidth',2,...
    'MarkerSize',10)
legend({'Empirical','Theoretical'},'FontSize',24,'Location','southeast')
grid on;
xlabel('\lambda', 'FontSize',24)
ylabel('MSE', 'FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)

figure(2)
hold on
plot(lambda_vec,varerror,...
    'LineWidth',2,...
    'MarkerSize',10)
grid on;
xlabel('\lambda', 'FontSize',24)
ylabel('Variance of Squared Error', 'FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)

figure(3)
hold on
plot(lambda_vec,savtwo,lambda_vec,savM,'--',lambda_vec,savone,':',...
    'LineWidth',2,...
    'MarkerSize',10)
legend({'Asymptotic Error','Effective Sparsity','Error L_1 norm'},'FontSize',24,'Location','southeast')
grid on;
xlabel('\lambda', 'FontSize',24)
ylabel('MSE and Sparsity', 'FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)
