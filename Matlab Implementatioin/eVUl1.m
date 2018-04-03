clc
clear
pr.sigma=0.01;
pr.epsilon =1e-6;
% A=[2,0;
%     0,1];
% b=[1;0.5];
%addpath .
[b,A]=libsvmread('../Datasets/bodyfat');
n=size(A,2);
mu=4;
x=0.1*ones(n,1);
M=A'*A;
U=zeros(n,n);
tau=0.1*norm(A'*b,Inf);
Dq = @(x) A'*(A*x-b);
fq=@(x) 0.5*(x'*A'-b')*(A*x-b);
f=@(x) fq(x)+tau*norm(x,1);
dqx=Dq(x);
p=zeros(n,1);
j=0;
while (1)
    alpha=tau/mu;
    w=x-1/mu *dqx;
    p(w>alpha)=w(w>alpha)-alpha;
    p(w<-alpha)=w(w<-alpha)+alpha;
    eps=f(x)-f(p)-mu*(p-x)'*(p-x);
    dqp=Dq(p);
    g=mu*(x-p)+dqp-dqx;
    for i=1:n
        if abs(x(i)) > eps / 2
            j=j+1;
            U(i,j)=1;
        end
    end
    if (j <1)
        
        xplus = p;
        
    else
        Q=U(:,j)'*M*U(:,j);
        ud=linsolve(Q,-U(:,j)'*g);
        xplus=p+U(:,j)*ud;
    end
    if f(x)-f(xplus)>=pr.sigma *(dqx'*(x-p)+tau*(norm(x,1)-norm(p,1)))
        x=xplus;
        dqx=Dq(x);
    else
        mu=mu*2;
    end
    j=0;
    U=zeros(n,n);
    p=zeros(n,1);
    if norm(g)<pr.epsilon
        break
    end
end
fv=min([f(x),f(xplus)]);
fprintf(1,'final value of the objective function = %6.3e\n',fv);