%To understand scaling, calculate the largest eigenvalues of the problem
%matrix A^T A.
addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
for j=1:l
    [b,A]=libsvmread(strcat(dataPath,nameFolds{j}));
    M=A'*A;
    maxEig=eigs(M,1);
   fprintf(nameFolds{j});
   fprintf(', max eigenvalue = %f\n',maxEig);
end