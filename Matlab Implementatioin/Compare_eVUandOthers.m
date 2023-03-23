%Completed comparision of eVU and VU with SpaRSA/m and gpsrBB!
clc;clear;
%parameters
pr.sigma=0.0001;pr.epsilon =1e-5;
pr.tmin=10*realmin;pr.tmax=realmax/10;
pr.eps_tol=1e-8;
pr.T=5;%max time
 pr.N=600;%max iteration
addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
fid=fopen('./OutputCompOthers1.txt','a');
fprintf(fid,'\nSave results of eVUl1 version 3.0 to compare with other 3 solvers ');
fprintf(fid,datestr(now));fprintf(fid,'\n');
fprintf(fid,'Parameters: pr.sigma=%f, pr.epsilon =%e, tau=0.1|ATb|_infty\n',pr.sigma,pr.epsilon);
fprintf(fid,'Starting point:  x=zeros(n,1)\n');
fprintf(fid,'Epsilon: epsilon=max([epsilon,0]);\n');
eVU_r=zeros(l,2);
VU_r=zeros(l,2);
admm=zeros(l,2);
for j=1:l
    [b,A]=libsvmread(strcat(dataPath,nameFolds{j}));
%     scaling
%     M=A'*A;
%     maxEig=eigs(M,1);
%     A=A/sqrt(maxEig);
%     b=b/sqrt(maxEig);
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
    %[fv,time]=eVUl1(x,M,f,Dq,tau,pr);
    addpath './v3/';
    [sol,history]=pqNVUlasso(x,M,f,Dq,tau,pr);
    fprintf(fid,'fv=%e, cpu=%f\n',sol.fv,sol.time+  telapsed0);
    eVU_r(j,:)=[sol.fv, sol.time + telapsed0];    
    addpath './V2.0/';
  [fv,time]=VUl1(x,M,f,Dq,tau,pr);
    fprintf(fid,'fv=%e, cpu=%f\n',fv,time+  telapsed0);
    VU_r(j,:)=[fv, time + telapsed0];
     addpath './admm/';
    [z, historya,time] = lasso(A, b, tau, 1.0, 1.0);
    admm(j,:)=[f(z),time];
end
save('results_compOther1.mat','eVU_r','VU_r');
difference= VU_r- eVU_r;
results1=sum(difference>0);
results2=sum(difference<0);
fprintf(fid,'Number of problems that eVU returns smaller value=%d, less CPU time=%d\n',results1(1),results1(2));
fprintf(fid,'Number of problems that eVU returns bigger value=%d, more CPU time=%d\n',results2(1),results2(2));
fprintf('fv(eVU)=fv(VU) for %d problems.',l-results1(1)-results2(1));
addpath 'C:\Users\Leo\Google Drive (shuai0liu@gmail.com)\0My writings\prox Newton\refs\SpaRSA_2.0';
load('otherResults.mat');
fv=[eVU_r(:,1), VU_r(:,1), resSpaRSA(:,1),resSpaRSA_m(:,1),resgpsrBB(:,1),admm(:,1)];
ha=perfprof(fv,1.05);
legend(ha,{'\(\varepsilon\)-\(\mathcal{VU}\)','\(\mathcal{VU}\)','SpaRSA','SpaRSAm','gpsrBB','admm'},'Location','southeast','Interpreter','latex');
%legend(ha,{'\(\varepsilon\)-\(\mathcal{VU}\)','\(\mathcal{VU}\)','SpaRSA','SpaRSAm','l1ls','gpsrBB'},'Location','southeast','Interpreter','latex');
figure;
cpu=[eVU_r(:,2), VU_r(:,2),resSpaRSA(:,2),resSpaRSA_m(:,2),resgpsrBB(:,2),admm(:,2)];
%cpu=[eVU_r(:,2), VU_r(:,2),resSpaRSA(:,2),resSpaRSA_m(:,2),resl1_ls(:,2),resgpsrBB(:,2)];
hb=perfprof(cpu,5);
legend(hb,{'\(\varepsilon\)-VU','\(\mathcal{VU}\)','SpaRSA','SpaRSAm','gpsrBB','admm'},'Location','southeast','Interpreter','latex');
%legend(hb,{'\(\varepsilon\)-VU','\(\mathcal{VU}\)','SpaRSA','SpaRSAm','l1ls','gpsrBB'},'Location','southeast','Interpreter','latex');
% ha=perfprof(fv,1.2);
% legend(ha,{'\(\varepsilon\)-\(\mathcal{VU}\)','\(\mathcal{VU}\)','SpaRSA','SpaRSAm','l1_ls','gpsrBB'},'Location','southeast','Interpreter','latex');
% fv=[eVU_r(:,1), resgpsrBB(:,1)];
% ha=perfprof(fv,2);
% legend(ha,{'\(\varepsilon\)-\(\mathcal{VU}\)','resgpsrBB'},'Location','southeast','Interpreter','latex');






