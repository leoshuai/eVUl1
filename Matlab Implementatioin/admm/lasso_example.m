% L1-regularized least-squares example

%% Generate problem data

% randn('seed', 0);
% rand('seed',0);
% 
% m = 1500;       % number of examples
% n = 5000;       % number of features
% p = 0.7;      % sparsity density  
% 
% x0 = sprandn(n,1,p);
% A = randn(m,n);
% A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
% b = A*x0 + sqrt(0.001)*randn(m,1);
load('../dataA.mat')
load('../datab.mat')
lambda_max = norm( A'*b, 'inf' );
lambda = 0.1*lambda_max; 
fq=@(x) 0.5*(x'*M*x-2*x'*atb+bb);
f=@(x) fq(x)+tau*norm(x,1);
%% Solve problem

[x, history] = lasso(A, b, lambda, 1.0, 1.0);

%% Reporting
K = length(history.objval);                                                                                                        
fbest=f(x);%1.244530977556365e+06
Y=(history.objval-fbest)/fbest;
h = figure;
plot(1:K, Y, 'k', 'MarkerSize', 10, 'LineWidth', 2); 
ylabel('f(x^k) -f(best)/f(best)'); xlabel('iter (k)');

% g = figure;
% subplot(2,1,1);                                                                                                                    
% semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
%     1:K, history.eps_pri, 'k--',  'LineWidth', 2); 
% ylabel('||r||_2'); 
% 
% subplot(2,1,2);                                                                                                                    
% semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
%     1:K, history.eps_dual, 'k--', 'LineWidth', 2);   
% ylabel('||s||_2'); xlabel('iter (k)'); 
