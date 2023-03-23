function [sol,history]=pqNVUlasso(x,M,f,Dq,tau,pr)
% V3.0
% 1. Return a solution and history
% 2.

tiv=tic;
n=length(x);
U = spalloc(n,n,n);
dqx=Dq(x);
p=zeros(n,1);
j=0;l=0;
k=0;%The iteration index
fx=f(x);
history.objval(k+1)=fx;
six=sign(x);
dfx=dqx+tau*six;
mu=0.5*(dfx'*dfx)/max(1,abs(fx));
%mu=max(pr.tmin,min(pr.tmax,mu));
%stop=false;
%tol_e=sqrt(eps);
while (true)
    alpha=tau/mu;
    w=x-(1/mu) *dqx;
    p(w>alpha)=w(w>alpha)-alpha;
    p(w<-alpha)=w(w<-alpha)+alpha;
    pmx=x-p;
    mx_mp=dqx'*pmx+tau*(norm(x,1)-norm(p,1));
    epsilon=mx_mp-mu*(pmx')*pmx;
    epsilon=max([epsilon,0]);
    %epsilon=max([epsilon,tol_e]);% This choice is not good!
    dqp=Dq(p);
    g=mu*pmx+dqp-dqx;
    nog=norm(g);
        if (nog<pr.epsilon)
            if  fx>fxpl
                sol.x=xplus;
                sol.fv=fxpl;
            else
                sol.x=x;
                sol.fv=fx;
            end
            sol.time=toc(tiv);
            fprintf('eVUl1 is terminated successfully. f_val=: %6.3e, time=%f\n',sol.fv,sol.time);
            return
        end
    for i=1:n
        if abs(x(i)) > epsilon / 2
            % if abs(x(i)) > 0
            j=j+1;
            U(i,j)=1;
        end
    end
    
    if (j <1)
        xplus = p;
    else
        Q=U(:,1:j)'*M*U(:,1:j);
        ud=mldivide(Q,-U(:,1:j)'*g);
        xplus=p+U(:,1:j)*ud;
    end
    fxpl=f(xplus);
    if fx-fxpl>=pr.sigma *mx_mp
        %Update mu
        dfx=dqx+tau*six;
        dxp=Dq(xplus);
        if j<1
            deltas=g-dfx;
        else
            deltas=dxp-dqx+tau*(sign(xplus)-six);
        end
        t1=deltas'*deltas;
        mu=t1*mu/(t1+mu*deltas'*(xplus-x));
        %mu=max(pr.tmin,min(pr.tmax,mu));
        x=xplus;
        dqx=dxp;
        k=k+1;
        fx=fxpl;
        history.objval(k+1)=fx;
    else
        % serious=false;
        mu=mu*2;
        % mu=min(mu*2,pr.tmax);
    end
    %       fprintf('k=%d, serious=%d, f(x)=%.4e, f(x+)=%.4e, epsilon=%.4e, norm(g)=%.4e, mu=%.3e,dim Uk=%d\n',...
    %        k,serious,fx,fxpl,epsilon,nog,mu,j);
    j=0;
    U=spalloc(n,n,n);
    p=zeros(n,1);
    l=l+1;
    temp=toc(tiv);
    %         switch pr.stopCri
    %             case 0
    %                 stop=(nog<pr.epsilon);
    %             case 1
    %                 stop= (min(fx,fxpl)<pr.epsilon );
    %             case 2
    %                 stop =(nog<pr.epsilon  && epsilon<pr.eps_tol);
    %         end
    %if (temp>pr.T || l>pr.N)
    if (temp>pr.T)
        %Considered fail
        if  fx>fxpl
            sol.x=xplus;
            sol.fv=fxpl;
        else
            sol.x=x;
            sol.fv=fx;
        end
        sol.time=toc(tiv);
        fprintf('f_val=: %6.3e, time=%f\n',sol.fv,sol.time);
        fprintf('l=: %d, mu=%e\n',l,mu);
        return
    end
end
end