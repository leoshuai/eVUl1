%V1 vs V2 for eVU 
% To be continued
clc;clear;
%parameters
pr.sigma=0.0001;pr.epsilon =1e-4;
pr.tmin=10*realmin;pr.tmax=realmax/10;
pr.eps_tol=1e-8;mu=4;pr.stopCri=0;
pr.T=2;%max time
 pr.N=300;%max iteration
pr.T=2;%max time
 pr.N=300;%max iteration
addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
fid=fopen('./Out_v1.txt','w');
fprintf(fid,'This is comparision of eVUl1 version 2.0 and v1. ');
fprintf(fid,datestr(now));fprintf(fid,'\n');
fprintf(fid,'Parameters: pr.sigma=%f, pr.epsilon =%e, mu=4, tau=0.1|ATb|_infty\n',pr.sigma,pr.epsilon);
fprintf(fid,'Starting point:  x=zeros(n,1)\n');
fprintf(fid,'Epsilon: epsilon=max([epsilon,0]) for eVU ;\n');
% eVU_r=zeros(l,2);
% VU_r=zeros(l,2);
addpath './V1.0/';
 load('eVU_VU.mat');
 VU_r=eVU_r;%VU_r contains eVU v2.0
 eVU_r=zeros(l,2);%eVU_r is v1.0
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
    [fv,time]=eVUl1(x,M,f,Dq,tau,pr,mu);
    fprintf(fid,'fv=%e, cpu=%f\n',fv,time+  telapsed0);
    eVU_r(j,:)=[fv, time + telapsed0];    
%   [fv,time]=VUl1(x,M,f,Dq,tau,pr);
%     fprintf(fid,'fv=%e, cpu=%f\n',fv,time+  telapsed0);
%     VU_r(j,:)=[fv, time + telapsed0];
end
save('eVU_V1.mat','eVU_r');%eVU_r is v1.0
difference= eVU_r- VU_r;
results2=sum(difference>0);
results1=sum(difference<0);
fprintf(fid,'Number of problems that eVU returns smaller value=%d, less CPU time=%d\n',results1(1),results1(2));
fprintf(fid,'Number of problems that eVU returns bigger value=%d, more CPU time=%d\n',results2(1),results2(2));
fprintf('fv(eVU)=fv(VU) for %d problems.',l-results1(1)-results2(1));

fv=[eVU_r(:,1), VU_r(:,1)];
ha=perfprof(fv,2);

legend(ha,{'\(\varepsilon\)-\(\mathcal{VU}\)','\(\mathcal{VU}\)'},'Location','southeast','Interpreter','latex');
figure;
cpu=[eVU_r(:,2), VU_r(:,2)];
hb=perfprof(cpu,2);
legend(hb,{'\(\varepsilon\)-VU','\(\mathcal{VU}\)'},'Location','southeast','Interpreter','latex');








