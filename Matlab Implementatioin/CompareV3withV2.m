load('r_compOtherV3.mat');
ver3res=VU_r;
load('results_compOther.mat');
ver2res=VU_r;
fv=[ver3res(:,1), ver2res(:,1)];
ha=perfprof(fv,2);
legend(ha,{'\(v3\)','\(v2\)'},'Location','southeast','Interpreter','latex');
figure;
cpu=[ver3res(:,2), ver2res(:,2)];
ha=perfprof(cpu,2);
legend(ha,{'\(v3\)','\(v2\)'},'Location','southeast','Interpreter','latex');


