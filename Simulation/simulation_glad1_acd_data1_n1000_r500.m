clc
clear all
load("simulation_glgd1_data1_n1000_r2000.mat")
randx = randx(2:1001,1:2000);
n=1000;
repeat = 2000;
alpha = 0.0793;
beta = 0.3046;
gamma = 0.6179;
phi01 = -0.8240;
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
    
    [parameters1, stderrors1, LLF1, mux1,~,~,~,pvalues1] = ACD_exponential(randxx',[],[]) ;%-786
    [parameters2, stderrors2, LLF2, mux2, ~,~,~,pvalues2] = ACD_weibull(randxx',[],[]);%-339
    [parameters3, stderrors3, LLF3, mux3, hx3,sx3,kx3,pvalues3] = ACD_generalized_gamma(randxx',[],[]);%-321
    [parameters4, stderrors4, LLF4, mux4, hx4,sx4,kx4,pvalues4] = ACD_burr(randxx',[],[]);% -316
    params1(re,:)=parameters1';
    params2(re,:)=parameters2';
    params3(re,:)=parameters3';
    params4(re,:)=parameters4';
    %%
    nu_len = 1;
    alpha0 = alpha;
    beta0 = beta;
    gamma0 = gamma;
    phi010 = phi01;
    x0 = [alpha0;beta0;gamma0;phi010];  
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(1,1)];
    upperbounds  = [1;1;1;10e8*ones(1,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)constant_loglikelihood1(params,randxx'),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=constant_loglikelihood1(xx,randxx');
    
    final_llf(1,nu_len)=llf;
    
    for i=2:10
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)constant_loglikelihood1(params,randxx'),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux5]=constant_loglikelihood1(xx,randxx);
        
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters5, stderrors5, ~, ~, ~, ~,~,~] = ACD_constant_nnsnp1(randxx',xx,[]) ;
    params5(re,:)=parameters5;
    sd5(re,:)=stderrors5;
end
save("simulation_glg1_acd_data1_n1000_r2000.mat")

    