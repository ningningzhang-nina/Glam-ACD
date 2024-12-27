clc
clear all
load("simulation_glgd1_data1_n1000_r2000.mat")
randx = randx(2:1001,1:1000);
n=1000;
repeat = 1000;
alpha = -0.2793;
beta = 0.2624;
gamma = 0.8701;
phi01 = -0.8240;
nu0=1;
for re=1:repeat
    psi = zeros(1,n);
    psi(1) = mean(randx(:,re));
    randxx = zeros(1,n);
    randxx(1)=randx(1,re)/exp(psi(1));
    for i = 2:n
        psi(i) = alpha + beta*randx(i-1,re) + gamma*psi(i-1);
        randxx(i) = exp(psi(i))*randx(i,re);
    end
    
    
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
        'objective',@(params)logconstant_loglikelihood1(params,randxx',2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood1(xx,randxx',2);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:10
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood1(params,randxx',2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux5]=logconstant_loglikelihood1(xx,randxx',2);
        
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters13, stderrors13, LLF13, mux13, hx13, pvalues13,AIC13,BIC13] = ACD_logconstant_nnsnp1(randxx',xx,2,[]) ;
    params5(re,:)=parameters13;
    sd5(re,:)=stderrors13;
end
save("simulation_glgd1_logacd2_data1_n1000_r1000.mat")

    