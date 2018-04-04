addpath C:\Users\Leo\Source\eVUl1;%path of libsvmread
dataPath='../Datasets/';% path of input datasets
d = dir(dataPath);nameFolds = {d(:).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
l=length(nameFolds);
% noisy observations
sigma = 0.0001;quiet=1;
for j=1:l
    [y,R]=libsvmread(strcat(dataPath,nameFolds{j}));
    % n is the original signal length
    n=size(R,2);
    % k is number of observations to make
    k =size(R,1);
    % number of  spikes to put down
    n_spikes = floor(.2*k);
    % n_spikes = 160;
    % set the handles for multiplying by R and R'
    hR = @(x) R*x;
    hRt = @(x) R'*x;
    % regularization parameter
    tau_max = max(abs(R'*y));
    tau = 0.1 * tau_max;
    fprintf(1,'\n\n-------------------------------------------------\n')
    fprintf(1,'-------------------------------------------------\n')
    fprintf(1,'Problem: n = %g,  k = %g, number of spikes = %g\n',n,k,n_spikes)
    fprintf(1,'Parameters: sigma = %g, tau = %g, debiasing = %g\n',sigma,tau,0)
    fprintf(1,'All BB algorithms initialized with zeros\n')
    fprintf(1,'-------------------------------------------------\n')
    
    
    [x_l1_ls,status,history] = l1_ls(R,y,2*tau,1.e-2,quiet);
    t_l1_ls = history(7,end);
    % for the reason explained above the objective
    % function of l1_ls is twice that of GPSR, so we
    % need the following correction.
    obj_l1_ls = history(2,end)/2.0;
    [dummy, iterations_l1_ls] = size(history);
    % now call fpc
    opts=fpc_opts([]);
    opts.mxitr=2000;
    opts.xtol=1.e-3;
    t0 = cputime;
    Out = fpc(n,R,y,1/tau,[],opts,[]);
    t_fpc = cputime - t0;
    x_fpc = Out.x;
    % need to scale function value:
    obj_fpc = tau*Out.f(end);
    iterations_fpc = Out.itr;
    tolA = obj_fpc;
    stopCri = 4;
    obj_BB_mono=[]; times_BB_mono=[];
    [x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
        times_BB_mono,debias_start_BB_mono,mse_BB_mono]= ...
        GPSR_BB(y,hR,tau,...
        'Debias',0,...
        'AT',hRt,...
        'Monotone',1,...
        'Initialization',0,...
        'StopCriterion',stopCri,...
        'ToleranceA',tolA,...
        'Verbose',0,...
        'Continuation',0);
    t_BB_mono = times_BB_mono(end);
    
    obj_SpaRSA_m=[]; times_SpaRSA_m=[];
    [x_SpaRSA_m,x_debias_SpaRSA_m,obj_SpaRSA_m,...
        times_SpaRSA_m,debias_start_SpaRSA_m,mse_SpaRSA_m,taus]= ...
        SpaRSA(y,hR,tau,...
        'Monotone',1,...
        'Debias',0,...
        'AT',hRt,...
        'Initialization',0,...
        'StopCriterion',stopCri,...
        'ToleranceA',tolA,...
        'MaxiterA',10000,...
        'Verbose',0,...
        'Continuation',0);
    t_SpaRSA_m = times_SpaRSA_m(end);
    
    obj_SpaRSA_s=[]; times_SpaRSA_s=[];
    [x_SpaRSA_s,x_debias_SpaRSA_s,obj_SpaRSA_s,...
        times_SpaRSA_s,debias_start_SpaRSA_s,mse_SpaRSA_s,taus]= ...
        SpaRSA(y,hR,tau,...
        'Monotone',0,...
        'Safeguard',1,...
        'M',5,...
        'Debias',0,...
        'AT',hRt,...
        'Initialization',0,...
        'StopCriterion',stopCri,...
        'ToleranceA',tolA,...
        'MaxiterA',10000,...
        'Verbose',0,...
        'Continuation',0);
    t_SpaRSA_s = times_SpaRSA_s(end);
    
    fprintf(1,'\n ------ WITHOUT CONTINUATION ------------------------ \n')
    
    fprintf(1,'SpaRSA-monotone; cpu: %6.2f secs (%d iterations)\n',...
        t_SpaRSA_m,length(obj_SpaRSA_m))
    fprintf(1,'final value of the objective function = %6.3e\n',...
        obj_SpaRSA_m(end));
     fid=fopen('./SpaRSA-monotone-output.txt','a');
     fprintf(fid,nameFolds{j});
     fprintf(fid,', fv=%e\n',obj_SpaRSA_m(end));
     
    fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_SpaRSA_m~=0))
    fprintf(1,'SpaRSA; cpu: %6.2f secs (%d iterations)\n',...
        t_SpaRSA_s,length(obj_SpaRSA_s))
    fprintf(1,'final value of the objective function = %6.3e\n',...
        obj_SpaRSA_s(end))
    fid=fopen('./SpaRSA-output.txt','a');
     fprintf(fid,nameFolds{j});
     fprintf(fid,', fv=%e\n',obj_SpaRSA_s(end));
    fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_SpaRSA_s~=0))
    
    fprintf(1,'\nGPSR-BB-monotone; cpu: %6.2f secs (%d iterations)\n',...
        t_BB_mono,length(obj_BB_mono))
    fprintf(1,'final value of the objective function = %6.3e\n',...
        obj_BB_mono(end))
    fid=fopen('./GPSR-BB-monotone-output.txt','a');
     fprintf(fid,nameFolds{j});
     fprintf(fid,', fv=%e\n',obj_BB_mono(end));
    fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_BB_mono~=0))
    
    fprintf(1,'\nl1_ls; cpu: %6.2f secs (%d iterations)\n',...
        t_l1_ls,iterations_l1_ls)
    fprintf(1,'final value of the objective function = %6.3e\n',...
        obj_l1_ls)
     fid=fopen('./l1_ls-output.txt','a');
     fprintf(fid,nameFolds{j});
     fprintf(fid,', fv=%e\n',obj_l1_ls);    
    fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_l1_ls~=0))
    
    fprintf(1,'\nFPC; cpu: %6.2f secs (%d iterations)\n',...
        t_fpc,iterations_fpc)
    fprintf(1,'final value of the objective function = %6.3e\n',...
        obj_fpc)
     fid=fopen('./FPC-output.txt','a');
     fprintf(fid,nameFolds{j});
     fprintf(fid,', fv=%e\n',obj_fpc);    
    fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_fpc~=0))
    fprintf(1,' ---------------------------------------------------- \n')
end




% Need to call l1_ls with tau/2 because it assumes
% a different objective:
%  || y - R*x||_2^2 + tau ||x||_1
% instead of the one assumed by GPSR which is
%  (1/2)*|| y - R*x||_2^2 + tau ||x||_1
%
