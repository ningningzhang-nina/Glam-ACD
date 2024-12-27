clc
clear all
load("simulation_glgd2_data2_n500_r1000.mat")
randx = randx(2:501,1:1000);
n=500;
repeat = 1000;
alpha = 0.1539;
beta = 0.2104;
gamma = 0.6404;
phi01 = -1.0188;
phi02=0.6696;
nu0=1;
for re=1:repeat
    psi = zeros(1,n);
    psi(1) = alpha/(1-beta-gamma);
    randxx = zeros(1,n);
    randxx(1)=alpha/(1-beta-gamma)*randx(1,re);
    for i = 2:n
        psi(i) = alpha + beta*randxx(i-1) + gamma*psi(i-1);
        randxx(i) = psi(i)*randx(i,re);
    end
    
    %%
    nu_len = 2;
    alpha0 = alpha;
    beta0 = beta;
    gamma0 = gamma;
    phi010 = phi01;
    phi020 = phi02;
    x0 = [alpha0;beta0;gamma0;phi010;phi020];  
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(2,1)];
    upperbounds  = [1;1;1;10e8*ones(2,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)constant_loglikelihood2(params,randxx'),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=constant_loglikelihood2(xx,randxx');
    
    final_llf(1,nu_len)=llf;
    
    for i=2:10
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)constant_loglikelihood2(params,randxx'),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux6]=constant_loglikelihood2(xx,randxx');
        xx
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters5, stderrors5, LLF5, mux5, hx5, sx5, kx5, pvalues5,AIC5,BIC5] = ACD_constant_nnsnp(randxx',xx,[]) ;

    params5(re,:)=parameters5';
    sd5(re,:)=stderrors5';
end
save("simulation_glgd2_acd_data2_n500_r1000.mat")