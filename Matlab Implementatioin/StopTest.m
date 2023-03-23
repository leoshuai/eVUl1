%Which stop criterion is good? 2 or 0?
%Conclusion: 0 is better!
clc;clear;
%parameters
pr.sigma=0.0001;pr.epsilon =1e-4;
pr.stopCri=2;pr.tmin=10*realmin;pr.tmax=realmax/10;
pr.eps_tol=1e-8;
pr.T=2;%max time
 pr.N=300;%max iteration
addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
%fi0=fopen('./Warnings.txt','a');
fid=fopen('./Output1.txt','w');
fprintf(fid,'This is comparision of stop 2 vs 0 for eVUl1 version 2.0 ');
%fprintf(fid,'Save results of eVUl1 version 2.0 and VUl1 Version 2.0 after removing exchanging 6 probs because VUl1 is stuck there ');
%fprintf(fid,'This is comparision of eVUl1 version 2.0 and VUl1 v2.0. ');
fprintf(fid,datestr(now));fprintf(fid,'\n');
fprintf(fid,'Parameters: pr.sigma=%f, pr.epsilon =%e, tau=0.1|ATb|_infty\n',pr.sigma,pr.epsilon);
fprintf(fid,'Starting point:  x=zeros(n,1)\n');
%fprintf(fid,'Epsilon: epsilon=max([epsilon,sqrt(eps)])  for eVU and no epsilon for VU;\n');
%fprintf(fid,'Epsilon: epsilon=max([epsilon,0]) for eVU and no epsilon for VU;\n');
fprintf(fid,'Epsilon: epsilon=max([epsilon,0]);\n');
fprintf(fid,'Stopping: 2 and then 0 \n');
eVU_re2=zeros(l,2);
eVU_re0=zeros(l,2);
addpath './V2.0/';
for j=1:l
    [b,A]=libsvmread(strcat(dataPath,nameFolds{j}));
    n=size(A,2);x=zeros(n,1);
    tstart = tic;
    M=A'*A;
    bb=b'*b;
    atb=A'*b;
    telapsed0 = toc(tstart);
%     To check if M is postive definite
%     r= rank(full(M));
%         if r~=n
%             fprintf(fi0,nameFolds{j});
%             fprintf(fi0, ': A^T A is not positive definite!\n');
%         end
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
    eVU_re2(j,:)=[fv, time + telapsed0];
end
%Now test VU with Stop=0.
pr.stopCri=0;
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
    eVU_re0(j,:)=[fv, time+  telapsed0];
end
%save('results.mat','eVU_re2','VU_re2');
difference= eVU_re2- eVU_re0;
results2=sum(difference>0);
results1=sum(difference<0);
fprintf(fid,'Number of problems that eVU2 returns smaller value=%d, less CPU time=%d\n',results1(1),results1(2));
fprintf(fid,'Number of problems that eVU2 returns bigger value=%d, more CPU time=%d\n',results2(1),results2(2));
fprintf('fv(eVU2)=fv(eVU0) for %d problems.',l-results1(1)-results2(1));
fv=[eVU_re2(:,1), eVU_re0(:,1)];
ha=perfprof(fv,2);
legend(ha,{'\(eVU2\)','\(eVU0\)'},'Location','southeast','Interpreter','latex');
figure;
cpu=[eVU_re2(:,2), eVU_re0(:,2)];
ha=perfprof(cpu,2);
legend(ha,{'\(eVU2\)','\(eVU0\)'},'Location','southeast','Interpreter','latex');
