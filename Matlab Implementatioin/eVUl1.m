function fv=eVUl1(x,M,f,Dq,tau,pr,mu)
xkm=x;
n=length(x);
U=zeros(n,n);
dqx=Dq(x);
p=zeros(n,1);
j=0;
k=0;%The iteration index
while (1)
    alpha=tau/mu;
    w=x-1/mu *dqx;
    p(w>alpha)=w(w>alpha)-alpha;
    p(w<-alpha)=w(w<-alpha)+alpha;
    mx_mp=dqx'*(x-p)+tau*(norm(x,1)-norm(p,1));
    eps=mx_mp-mu*(p-x)'*(p-x);
    fprintf('\nIter. %d: epsilon =%e; ',k,eps);
    if eps<0
        fprintf('\nWarning: epsilon <0, something is wrong!\n');
    end
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
       % fprintf('Warning: Uk  is vacuous \n');
    else
        Q=U(:,j)'*M*U(:,j);
        ud=linsolve(Q,-U(:,j)'*g);
        xplus=p+U(:,j)*ud;
        fprintf(' dim Uk=%d \n',j);
    end
    if f(x)-f(xplus)>=pr.sigma *mx_mp
        %Update mu
        t=max([0.04, 0.1*mu]);
        gxkm=Dq(xkm)+tau*Dh(xkm);
        t1=Dq(x)+tau*Dh(x)-gxkm;
        mun=t1'*t1;
        mun=mun/(t1'*(x-xkm));
        t=max([mun,t]);
        mu=min([10*mu,t]);
        xkm=x;
        x=xplus;
        dqx=Dq(x);
        k=k+1;
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
fprintf(1,'eVUl1 is terminated successfully with final value of the objective function: %6.3e\n',fv);
end