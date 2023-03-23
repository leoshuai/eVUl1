%Which setting for epsilon?epsilon=max([epsilon,0]) or
%max([epsilon,sqrt(eps)])?
%Conclusion: the first  is better!
%ALso did test with stop=2 and max([epsilon,sqrt(eps)]). not good.
clc;clear;
%parameters
pr.sigma=0.0001;pr.epsilon =1e-4;
%pr.stopCri=2;
pr.tmin=10*realmin;pr.tmax=realmax/10;
pr.eps_tol=1e-8;
pr.T=2;%max time
 pr.N=300;%max iteration
addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
fid=fopen('./Output2.txt','w');
fprintf(fid,'Which setting for epsilon?epsilon=max([epsilon,0]) or max([epsilon,sqrt(eps)])? for eVUl1 version 2.0. ');
fprintf(fid,'\nFirst with max([epsilon,sqrt(eps)]);the other set (with Stop=0 and max([epsilon,0]) )is already tested in Stoptest');
fprintf(fid,datestr(now));fprintf(fid,'\n');
fprintf(fid,'Parameters: pr.sigma=%f, pr.epsilon =%e, tau=0.1|ATb|_infty\n',pr.sigma,pr.epsilon);
fprintf(fid,'Starting point:  x=zeros(n,1)\n');
fprintf(fid,'Epsilon: epsilon=max([epsilon,sqrt(eps)]);\n');
eVU_s=zeros(l,2);
addpath './V2.0/';
for j=1:l
    [b,A]=libsvmread(strcat(dataPath,nameFolds{j}));
    n=size(A,2);x=zeros(n,1);
    tstart = tic;
    M=A'*A;
    bb=b'*b;
    atb=A'*b;
    telapsed0 = toc(tstart);
    tau=0.1*norm(atb,Inf);
    Dq = @(x) M*x-atb;
    fq=@(x) 0.5*(x'*M*x-2*x'*atb+bb);
    f=@(x) fq(x)+tau*norm(x,1);
    fprintf(fid,nameFolds{j});
    fprintf(fid,', n=%d\n ',n);
    fprintf(nameFolds{j});
    fprintf(', n=%d\n',n);
    [fv,time]=eVUl1(x,M,f,Dq,tau,pr);
    fprintf(fid,'fv=%e, cpu=%f\n',fv,time+  telapsed0);
    eVU_s(j,:)=[fv, time + telapsed0];
end
save('reEpsilonTest.mat','eVU_s');
%Now test eVU with Stop=0.
% eVU_re0=zeros(l,2);
% for j=1:l
%     [b,A]=libsvmread(strcat(dataPath,nameFolds{j}));
%     n=size(A,2);x=zeros(n,1);
%     tstart = tic;
%     M=A'*A;
%     bb=b'*b;
%     atb=A'*b;
%     telapsed0 = toc(tstart);
%     tau=0.1*norm(atb,Inf);
%     Dq = @(x) M*x-atb;
%     fq=@(x) 0.5*(x'*M*x-2*x'*atb+bb);
%     f=@(x) fq(x)+tau*norm(x,1);
%     fprintf(fid,nameFolds{j});
%     fprintf(fid,', n=%d\n ',n);
%     fprintf(nameFolds{j});
%     fprintf(', n=%d\n',n);
%     [fv,time]=eVUl1(x,M,f,Dq,tau,pr);
%     fprintf(fid,'fv=%e, cpu=%f\n',fv,time+  telapsed0);
%     eVU_re0(j,:)=[fv, time+  telapsed0];
% end
load('results.mat');
difference= eVU_s- eVU_re0;
results2=sum(difference>0);
results1=sum(difference<0);
fprintf(fid,'Number of problems that eVUs returns smaller value=%d, less CPU time=%d\n',results1(1),results1(2));
fprintf(fid,'Number of problems that eVUs returns bigger value=%d, more CPU time=%d\n',results2(1),results2(2));
fprintf('fv(eVUs)=fv(eVU0) for %d problems.',l-results1(1)-results2(1));
fv=[eVU_s(:,1), eVU_re0(:,1)];
ha=perfprof(fv,2);
legend(ha,{'\(eVUs\)','\(eVU0\)'},'Location','southeast','Interpreter','latex');
figure;
cpu=[eVU_s(:,2), eVU_re0(:,2)];
ha=perfprof(cpu,2);
legend(ha,{'\(eVUs\)','\(eVU0\)'},'Location','southeast','Interpreter','latex');
