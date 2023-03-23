% load('results.mat');
% load('otherResults.mat');
%A=[eVU_results(:,2), VU_results(:,2),resSpaRSA_m(:,2),resSpaRSA(:,2),resgpsrBB(:,2),resl1_ls(:,2),resfpc(:,2)];
%  cpu=[eVU_results(:,2), resSpaRSA(:,2)];
%  ha=perfprof(cpu,3);
% % legend(ha,'eVU','VU','SpaRSA_m','SpaRSA','GPSR','l1_ls','fpc');
%  legend(ha,'eVU','SpaRSA');
% F=[eVU_results(:,1), VU_results(:,1),resSpaRSA_m(:,1),resSpaRSA(:,1),resgpsrBB(:,1),resl1_ls(:,1),resfpc(:,1)];
% hf=perfprof(F);
% legend(hf,'eVU','VU','SpaRSA_m','SpaRSA','GPSR','l1_ls','fpc');
%  F=[eVU_results(:,1),resSpaRSA_m(:,1)];
%  hf=perfprof(F,1.2);
%  legend(hf,'eVU','SpaRSA_m');
%  difference= resSpaRSA- eVU_results;
% results1=sum(difference>0);
% results2=sum(difference<0);
load('results2.mat');
RStop2=eVU_results;
load('results.mat');
RStop0=eVU_results;
% cpu=[RStop2(:,2), RStop0(:,2)];
% ha=perfprof(cpu);
%  legend(ha,'RStop2','RStop0');
 % Camparing v1 and v2, we found v2 is better but eVU and VU still are
 % similar.
 
 
 fv=[RStop2(:,1), RStop0(:,1)];
ha=perfprof(fv);
 legend(ha,'RStop2','RStop0');
 