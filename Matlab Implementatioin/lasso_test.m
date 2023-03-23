% L1-regularized least-squares example
clc;clear;
%% Generate problem data

rng('default');
m = 2000;       % number of examples
n = 2000;       % number of features n>m
density = 0.8;      % sparsity density
np=10;% instances
raadmm=zeros(np,1);
rapqnvu=zeros(np,1);
i=0;
pr.sigma=0.0001;pr.epsilon =1e-6;
pr.tmin=10*realmin;pr.tmax=realmax/10;
pr.eps_tol=1e-8;
pr.T=2;%max time
pr.N=300;%max iteration
while (true)
    x0 = sprandn(n,1,density);
    A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    b = A*x0 + sqrt(0.001)*randn(m,1);
     %b = -20+100*rand(m,1);
    tstart = tic;
    M=A'*A;
    bb=b'*b;
    atb=A'*b;
    telapsed0 = toc(tstart);
    [~,p] = chol(M);
    %check if the U-Hessian is positive definite
    %because I cannot determine U(\bar x), I just check
    % if M is positive definite.
    tau=0.1*norm(atb,Inf);
    Dq = @(x) M*x-atb;
    fq=@(x) 0.5*(x'*M*x-2*x'*atb+bb);
    f=@(x) fq(x)+tau*norm(x,1);
    
    %% parameters
    x0 = zeros(n,1);
    %% Solve problem
    addpath './v3/';
    [sol,history]=pqNVUlasso(x0,M,f,Dq,tau,pr);
    addpath './admm/';
    [z, historya] = lasso(A, b, tau, 1.0, 1.0);
    fva=f(z);
    dif=(fva-sol.fv)
    if abs(dif)<1e-5
        i=i+1;
        disp('\nless than 1e-5\n')
        fbar=min(fva,sol.fv);
        raadmm(i)=-log10(max(10^(-16),(fva-fbar)/(1+abs(fbar))));
        rapqnvu(i)=-log10(max(10^(-16),(sol.fv-fbar)/(1+abs(fbar))));
    else
        disp('\ngreater than 1e-5\n')
    end
    if i==np
        break
    end
end

%% Reporting
save('./results/RA.mat','raadmm','rapqnvu');

% b = -20+100*rand(m,1);
% R = sprandsym(m,density,100*rand(m,1));
% [~,p] = chol(R);% verify R is positive definite;
%%GramSchmidt
% V=R(:,1:n);
%     A = zeros(m,n);
%     A(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
%     for i = 2:n
%       A(:,i) = V(:,i);
%       for j = 1:i-1
%         A(:,i) = A(:,i) - ( A(:,i)'*A(:,j) )/( A(:,j)'*A(:,j) )*A(:,j);
%       end
%       A(:,i) = A(:,i)/sqrt(A(:,i)'*A(:,i));
%     end
%     %A is approxiamtely orthonormal
% I=eye(n);
% dif2=I-A'*A;
% e=norm(dif2); %the error of orthonormality of A.
% A=R(:,1:n);
% save('dataA.mat','A');
% load('dataA.mat','A');
% load('datab.mat','b');
% In order to make sure 0 \in \ri \partial f(\bar x), we
% require A'A=I and exist an i such that the i-th component of
%c=(A'b/ tau) is in (-1,1)
% Thus M is alwayd positive definite
% c=atb/tau;%Check c has component in (-1,1).
% for i=1:m
%     if c(i)>-1 && c(i)<1
%         disp('0 \in \ri \partial f(\bar x)');
%         break
%     end
% end
% fid=fopen('./results/lassoTest.txt','a');
% fprintf(fid,datestr(now));fprintf(fid,'\n');
% addpath './V2.0/';
% [fv2,time2]=eVUl1(x0,M,f,Dq,tau,pr)
% K = length(history.objval);
% fbest=history.objval(end);%1.011798667697997e+06
% Y=(history.objval-fbest)/fbest;
% h = figure;
% plot(1:K, Y, 'k', 'MarkerSize', 10, 'LineWidth', 2);
% ylabel('(f(x^k)-f(best))/f(best)'); xlabel('iter (k)');
% fbest=1.244530977556365e+06;
% fbar=1.011798667697997e+06;