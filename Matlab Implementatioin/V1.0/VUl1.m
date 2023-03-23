function fv=VUl1(x,M,f,Dq,tau,pr,mu)
xkm=x;
n=length(x);
U=sparse(n,n);% Later use sparse
dqx=Dq(x);
p=zeros(n,1);
j=0;
k=0;%The iteration index
fx=f(x);
while (1)
    alpha=tau/mu;
    w=x-(1/mu) *dqx;
    p(w>alpha)=w(w>alpha)-alpha;
    p(w<-alpha)=w(w<-alpha)+alpha;
    mx_mp=dqx'*(x-p)+tau*(norm(x,1)-norm(p,1));
    %     mx_mp- 1537155252
    %     t=mu*(p-x)'*(p-x);
    %     t- 1537155252
    epsilon=mx_mp-mu*(p-x)'*(p-x);
    %     fprintf('epsilon =%f, ',epsilon);
    %     if epsilon<0
    %         fprintf('\nWarning: epsilon <0, something is wrong!\n');
    %     end
    dqp=Dq(p);
    g=mu*(x-p)+dqp-dqx;
    nog=norm(g);
    for i=1:n
        if abs(x(i)) > 0
            j=j+1;
            U(i,j)=1;
        end
    end
    if (j <1)
        xplus = p;
    else
        Q=U(:,j)'*full(M)*U(:,j);
        ud=linsolve(Q,-U(:,j)'*g);
        xplus=p+U(:,j)*ud;
    end
    fxpl=f(xplus);
    if fx-fxpl>=pr.sigma *mx_mp
        %Update mu
        serious=true;
        t=max([0.04, 0.1*mu]);
        gxkm=Dq(xkm)+tau*Dh(xkm);
        t1=dqx+tau*Dh(x)-gxkm;
        mun=t1'*t1;
        mun=mun/(t1'*(x-xkm));
        t=max([mun,t]);
        mu=min([10*mu,t]);
        xkm=x;
        x=xplus;
        fx=f(x);
        dqx=Dq(x);
        k=k+1;
    else
        serious=false;
        mu=mu*2;
    end
    fprintf('k=%d, serious=%d, f(x)=%.4e, f(x+)=%.4e,epsilon=%.4e, norm(g)=%.4e, mu=%.3e,dim Uk=%d\n',...
        k,serious,fx,fxpl,epsilon,nog,mu,j);
    j=0;
    U=zeros(n,n);
    p=zeros(n,1);
    if nog<pr.epsilon && epsilon<1e-8
        break
    end
end
fv=min([fx,fxpl]);
fprintf('VUl1 is terminated successfully. f_val=: %6.3e\n',fv);
end