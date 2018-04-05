clc;clear;
% A=[2,0,0;
%     0,1,0;
%     0,0,0.01];
% b=[1;0.5;0.1];
%parameters
digits(40)
pr.sigma=0.0001;pr.epsilon =1e-4;mu=4;
addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
for j=1:2
    
    [b,A]=libsvmread(strcat(dataPath,nameFolds{j}));
    n=size(A,2);x=zeros(n,1);M=A'*A;
    bb=b'*b;r= rank(full(M));
    % To check if M is postive definite
    if r~=n
        disp('Warning: A^T A is not positive definite!');
    end
    atb=A'*b;
    tau=0.1*norm(atb,Inf);
    Dq = @(x) M*x-atb;
    fq=@(x) 0.5*(x'*M*x-2*x'*atb+bb);
    f=@(x) fq(x)+tau*norm(x,1);
    fid=fopen('./output.txt','a');
     fprintf(fid,nameFolds{j});
     fv=eVUl1(x,M,f,Dq,tau,pr,mu);
     fprintf(fid,', fv=%e\n',fv);
end