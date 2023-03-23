function [fv,time]=eVUl1(x,M,f,Dq,tau,pr)
% V2.0
% Changed the update formula of mu
tiv=tic;
n=length(x);
U=sparse(n,n);% Later use sparse
dqx=Dq(x);
p=zeros(n,1);
j=0;l=0;
k=0;%The iteration index
fx=f(x);
dfx=dqx+tau*sign(x);
mu=0.5*(dfx'*dfx)/max(1,abs(fx));
mu=max(pr.tmin,min(pr.tmax,mu));
stop=false;
tol_e=sqrt(eps);
while (~stop)
    alpha=tau/mu;
    w=x-(1/mu) *dqx;
    p(w>alpha)=w(w>alpha)-alpha;
    p(w<-alpha)=w(w<-alpha)+alpha;
    mx_mp=dqx'*(x-p)+tau*(norm(x,1)-norm(p,1));
    epsilon=mx_mp-mu*(p-x)'*(p-x);
    %epsilon=max([epsilon,0]);
    epsilon=max([epsilon,tol_e]);
    %     fprintf('epsilon =%f, ',epsilon);
    %     if epsilon<0
    %         fprintf('\nWarning: epsilon <0, something is wrong!\n');
    %     end
    dqp=Dq(p);
    g=mu*(x-p)+dqp-dqx;
    nog=norm(g);
    for i=1:n
        if abs(x(i)) > epsilon / 2
            j=j+1;
            U(i,j)=1;
        end
    end
    
    if (j <1)
        xplus = p;
    else
        Q=U(:,j)'*M*U(:,j);
        %ud=linsolve(Q,-U(:,j)'*g);
        ud=mldivide(Q,-U(:,j)'*g);
        xplus=p+U(:,j)*ud;
    end
    fxpl=f(xplus);
    if fx-fxpl>=pr.sigma *mx_mp
        %Update mu
        % serious=true;
        
        dfx=dqx+tau*sign(x);
        if j<1
            deltas=g-dfx;
        else
            deltas=Dq(xplus)-dqx+tau*(sign(xplus)-sign(x));
        end
        t1=deltas'*deltas;
        mu=t1*mu/(t1+mu*deltas'*(xplus-x));
        mu=max(pr.tmin,min(pr.tmax,mu));
        x=xplus;
        fx=f(x);
        dqx=Dq(x);
        k=k+1;
    else
        % serious=false;
        mu=mu*2;
        % mu=min(mu*2,pr.tmax);
    end
    %       fprintf('k=%d, serious=%d, f(x)=%.4e, f(x+)=%.4e, epsilon=%.4e, norm(g)=%.4e, mu=%.3e,dim Uk=%d\n',...
    %        k,serious,fx,fxpl,epsilon,nog,mu,j);
    j=0;
    U=sparse(n,n);
    p=zeros(n,1);
    l=l+1;
    temp=toc(tic);
    stop=(nog<pr.epsilon);
    %         switch pr.stopCri
    %             case 0
    %                 stop=(nog<pr.epsilon);
    %             case 1
    %                 stop= (min(fx,fxpl)<pr.epsilon );
    %             case 2
    %                 stop =(nog<pr.epsilon  && epsilon<pr.eps_tol);
    %         end
    if (temp>pr.T || l>pr.N) && (stop==false)
        %Considered fail
        fv=min([fx,fxpl]);
        time=toc(tiv);
        fail=1;
        return
    end
end
fv=min([fx,fxpl]);
time=toc(tiv);
fprintf('eVUl1 is terminated successfully. f_val=: %6.3e, time=%f\n',fv,time);
end