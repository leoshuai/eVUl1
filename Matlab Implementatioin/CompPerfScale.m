clc
clear
% load('results_Scale.mat');
% addpath 'C:\Users\Leo\Google Drive (shuai0liu@gmail.com)\0My writings\Eps VU Algorithm\refs\SpaRSA_2.0';
% load('otherResultsScale.mat');
 cpu=[eVU_results(:,2), resSpaRSA(:,2)];
 ha=perfprof(cpu,2);
 legend(ha,'eVU','SpaRSA');
 
%  F=[eVU_results(:,1), resSpaRSA(:,1)];
%  hf=perfprof(F,1.2);
%  legend(hf,'eVU','SpaRSA');