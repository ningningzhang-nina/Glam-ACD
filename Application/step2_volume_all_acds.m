clc
clear all
%names ={'ABT','AMD','XOM','VZ','PM','PG','GOOGL','GE','DVN','DIS','CVX','ORCL','MRK','LOW','KO','JNJ'};
names = {'BUD','ADBE','AZN','LMT','SHW','SCCO'};
%load('volume_all_acd_results.mat')
% final_table = readtable('./new2024/BUD_volume_std.xlsx');
% %% step 1
% sta_durations=final_table.std_duration;
% inx=sta_durations(1:round(length(sta_durations)*2/3-1),1);
% outx=sta_durations(round(length(sta_durations)*2/3):end,1);
% outx1=sta_durations(round(length(sta_durations)*2/3-2):end-1,1);
for item =1:length(names)
    final_table = readtable(sprintf('./new2024/%s_volume_std.xlsx',names{item}));
    %% step 1
    sta_durations=final_table.std_duration;
    inx=sta_durations(1:round(length(sta_durations)*2/3-1),1);
    outx=sta_durations(round(length(sta_durations)*2/3):end,1);
    outx1=sta_durations(round(length(sta_durations)*2/3-2):end-1,1);
    %%
    [parameters1, stderrors1, LLF1, mux1,~,~,~,pvalues1] = ACD_exponential(inx,[],[]) ;%-786
    final_llf(item,1)=LLF1;
    params{item,1} = [parameters1,stderrors1,pvalues1];
    %%
    inz = insamplez(inx,mux1','eacd',parameters1);
    
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux1(end),'eacd',parameters1);
    inzs{item}(:,1) = inz;
    outzs{item}(:,1) = outz;
    
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(1,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %%
    [parameters2, stderrors2, LLF2, mux2, ~,~,~,pvalues2] = ACD_weibull(inx,[],[]);%-339
    final_llf(item,2)=LLF2;
    params{item,2} = [parameters2,stderrors2,pvalues2];
    %%
    inz = insamplez(inx,mux2','wacd',parameters2);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux2(end),'wacd',parameters2);
    inzs{item}(:,2) = inz;
    outzs{item}(:,2) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(2,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %% 
    [parameters3, stderrors3, LLF3, mux3, hx3,sx3,kx3,pvalues3] = ACD_generalized_gamma(inx,[],[]);%-321
    final_llf(item,3)=LLF3;
    params{item,3} = [parameters3,stderrors3,pvalues3];
    %%
    inz = insamplez(inx,mux3','ggacd',parameters3);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux3(end),'ggacd',parameters3);
    inzs{item}(:,3) = inz;
    outzs{item}(:,3) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(3,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %% 
    [parameters4, stderrors4, LLF4, mux4, hx4,sx4,kx4,pvalues4] = ACD_burr(inx,[],[]);% -316
    final_llf(item,4)=LLF4;
    params{item,4} = [parameters4,stderrors4,pvalues4];
    %%
    inz = insamplez(inx,mux4','bacd',parameters4);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux4(end),'bacd',parameters4);
    inzs{item}(:,4) = inz;
    outzs{item}(:,4) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(4,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %%
    nu_len = 1;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    x0 = [alpha;beta;gamma;phi01];  
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(1,1)];
    upperbounds  = [1;1;1;10e8*ones(1,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)constant_loglikelihood1(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=constant_loglikelihood1(xx,inx);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)constant_loglikelihood1(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux5]=constant_loglikelihood1(xx,inx);
        
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters5, stderrors5, LLF5, mux5, hx5, pvalues5,AIC5,BIC5] = ACD_constant_nnsnp1(inx,xx,[]);
    params{item,5} = [parameters5,stderrors5,pvalues5];
    final_llf(item,5)=llf;
    inz = insamplez(inx,mux5','led1',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux5(end),'led1',xx);
    inzs{item}(:,5) = inz;
    outzs{item}(:,5) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(5,:)=[inp,in2p,outp,out2p,inpp,outpp];
    
    nu_len = 2;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02];  
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(2,1)];
    upperbounds  = [1;1;1;10e8*ones(2,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)constant_loglikelihood2(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=constant_loglikelihood2(xx,inx);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)constant_loglikelihood2(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux6]=constant_loglikelihood2(xx,inx);
        xx
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters6, stderrors6, LLF6, mux6, hx6,sx6,kx6, pvalues6,AIC6,BIC6] = ACD_constant_nnsnp(inx,xx,[]) ;
    params{item,6} = [parameters6,stderrors6,pvalues6];
    final_llf(item,6)=llf;
    inz = insamplez(inx,mux6','led2',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux6(end),'led2',xx);
    inzs{item}(:,6) = inz;
    outzs{item}(:,6) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(6,:)=[inp,in2p,outp,out2p,inpp,outpp];
    nu_len = 3;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02;phi03]; 
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(3,1)];
    upperbounds  = [1;1;1;10e8*ones(3,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)constant_loglikelihood3(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=constant_loglikelihood3(xx,inx);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)constant_loglikelihood3(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux7]=constant_loglikelihood3(xx,inx);
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters7, stderrors7, LLF7, mux7, hx7,sx7,kx7, pvalues7,AIC7,BIC7] = ACD_constant_nnsnp3(inx,xx,[]) ;
    params{item,7} = [parameters7,stderrors7,pvalues7];
    final_llf(item,7)=llf;
    inz = insamplez(inx,mux7','led3',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux7(end),'led3',xx);
    inzs{item}(:,7) = inz;
    outzs{item}(:,7) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(7,:)=[inp,in2p,outp,out2p,inpp,outpp];
    nu_len = 4;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    phi04 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02;phi03;phi04];    
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(4,1)];
    upperbounds  = [1;1;1;10e8*ones(4,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)constant_loglikelihood4(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=constant_loglikelihood4(xx,inx);
    
    %final_llf(1,nu_len)=llf;
    if item==3
        whole_number=1000;
    else
        whole_number=100;
    end
    for i=2:whole_number
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)constant_loglikelihood4(params,inx),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux8]=constant_loglikelihood4(xx,inx);
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters8, stderrors8, LLF8, mux8, hx8,sx8,kx8, pvalues8,AIC8,BIC8] = ACD_constant_nnsnp4(inx,xx,[]) ;
    params{item,8} = [parameters8,stderrors8,pvalues8];
    final_llf(item,8)=llf;
    inz = insamplez(inx,mux8','led4',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux8(end),'led4',xx);
    inzs{item}(:,8) = inz;
    outzs{item}(:,8) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(8,:)=[inp,in2p,outp,out2p,inpp,outpp];


    %%
    [parameters1, stderrors1, LLF1, mux1,~,~,~,pvalues1] = ACD_logexponential(inx,[],1,[]) ;%-786
    final_llf(item,9)=LLF1;
    params{item,9} = [parameters1,stderrors1,pvalues1];
    %%
    inz = insamplez(inx,mux1','eacd',parameters1);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux1(end),'logeacd1',parameters1);
    inzs{item}(:,9) = inz;
    outzs{item}(:,9) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(9,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %%
    [parameters2, stderrors2, LLF2, mux2, ~,~,~,pvalues2] = ACD_logweibull(inx,[],1,[]);%-339
    final_llf(item,10)=LLF2;
    params{item,10} = [parameters2,stderrors2,pvalues2];
    %%
    inz = insamplez(inx,mux2','wacd',parameters2);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux2(end),'logwacd1',parameters2);
    inzs{item}(:,10) = inz;
    outzs{item}(:,10) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(10,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %% 
    [parameters3, stderrors3, LLF3, mux3, hx3,sx3,kx3,pvalues3] = ACD_loggeneralized_gamma(inx,[],1,[]);%-321
    final_llf(item,11)=LLF3;
    params{item,11} = [parameters3,stderrors3,pvalues3];
    %%
    inz = insamplez(inx,mux3','ggacd',parameters3);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux3(end),'logggacd1',parameters3);
    inzs{item}(:,11) = inz;
    outzs{item}(:,11) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(11,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %% 
    [parameters4, stderrors4, LLF4, mux4, hx4,sx4,kx4,pvalues4] = ACD_logburr(inx,[],1,[]);% -316
    final_llf(item,12)=LLF4;
    params{item,12} = [parameters4,stderrors4,pvalues4];
    %%
    inz = insamplez(inx,mux4','bacd',parameters4);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux4(end),'logbacd1',parameters4);
    inzs{item}(:,12) = inz;
    outzs{item}(:,12) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(12,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %%
    nu_len = 1;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    x0 = [alpha;beta;gamma;phi01];  
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(1,1)];
    upperbounds  = [1;1;1;10e8*ones(1,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood1(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood1(xx,inx,1);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood1(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux5]=logconstant_loglikelihood1(xx,inx,1);
        
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters13, stderrors13, LLF13, mux13, hx13, pvalues13,AIC13,BIC13] = ACD_logconstant_nnsnp1(inx,xx,1,[]) ;
    params{item,13} = [parameters13,stderrors13,pvalues13];
    final_llf(item,13)=llf;
    inz = insamplez(inx,mux5','led1',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux5(end),'logled11',xx);
    inzs{item}(:,13) = inz;
    outzs{item}(:,13) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(13,:)=[inp,in2p,outp,out2p,inpp,outpp];
    
    nu_len = 2;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02]; 
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(2,1)];
    upperbounds  = [1;1;1;10e8*ones(2,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood2(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood2(xx,inx,1);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood2(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux6]=logconstant_loglikelihood2(xx,inx,1);
        xx
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters14, stderrors14, LLF14, mux14, hx14,sx14,kx14, pvalues14,AIC14,BIC14] = ACD_logconstant_nnsnp(inx,xx,1,[]) ;
    params{item,14} = [parameters14,stderrors14,pvalues14];
    final_llf(item,14)=llf;
    inz = insamplez(inx,mux6','led2',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux6(end),'logled21',xx);
    inzs{item}(:,14) = inz;
    outzs{item}(:,14) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(14,:)=[inp,in2p,outp,out2p,inpp,outpp];
    nu_len = 3;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02;phi03]; 
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(3,1)];
    upperbounds  = [1;1;1;10e8*ones(3,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood3(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood3(xx,inx,1);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood3(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux7]=logconstant_loglikelihood3(xx,inx,1);
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters15, stderrors15, LLF15, mux15, hx15,sx15,kx15, pvalues15,AIC15,BIC15] = ACD_logconstant_nnsnp3(inx,xx,1,[]) ;
    params{item,15} = [parameters15,stderrors15,pvalues15];
    final_llf(item,15)=llf;
    inz = insamplez(inx,mux7','led3',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux7(end),'logled31',xx);
    inzs{item}(:,15) = inz;
    outzs{item}(:,15) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(15,:)=[inp,in2p,outp,out2p,inpp,outpp];
    nu_len = 4;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    phi04 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02;phi03;phi04];    
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(4,1)];
    upperbounds  = [1;1;1;10e8*ones(4,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood4(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood4(xx,inx,1);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood4(params,inx,1),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux8]=logconstant_loglikelihood4(xx,inx,1);
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters16, stderrors16, LLF16, mux16, hx16,sx16,kx16, pvalues16,AIC16,BIC16] = ACD_logconstant_nnsnp4(inx,xx,1,[]) ;
    params{item,16} = [parameters16,stderrors16,pvalues16];
    final_llf(item,16)=llf;
    inz = insamplez(inx,mux8','led4',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux8(end),'logled41',xx);
    inzs{item}(:,16) = inz;
    outzs{item}(:,16) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(16,:)=[inp,in2p,outp,out2p,inpp,outpp];

    %%
    [parameters1, stderrors1, LLF1, mux1,~,~,~,pvalues1] = ACD_logexponential(inx,[],2,[]) ;%-786
    final_llf(item,17)=LLF1;
    params{item,17} = [parameters1,stderrors1,pvalues1];
    %%
    inz = insamplez(inx,mux1','eacd',parameters1);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux1(end),'logeacd2',parameters1);
    inzs{item}(:,17) = inz;
    outzs{item}(:,17) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(17,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %%
    [parameters2, stderrors2, LLF2, mux2, ~,~,~,pvalues2] = ACD_logweibull(inx,[],2,[]);%-339
    final_llf(item,18)=LLF2;
    params{item,18} = [parameters2,stderrors2,pvalues2];
    %%
    inz = insamplez(inx,mux2','wacd',parameters2);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux2(end),'logwacd2',parameters2);
    inzs{item}(:,18) = inz;
    outzs{item}(:,18) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(18,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %% 
    [parameters3, stderrors3, LLF3, mux3, hx3,sx3,kx3,pvalues3] = ACD_loggeneralized_gamma(inx,[],2,[]);%-321
    final_llf(item,19)=LLF3;
    params{item,19} = [parameters3,stderrors3,pvalues3];
    %%
    inz = insamplez(inx,mux3','ggacd',parameters3);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux3(end),'logggacd2',parameters3);
    inzs{item}(:,19) = inz;
    outzs{item}(:,19) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(19,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %% 
    if item==4
        [parameters4, stderrors4, LLF4, mux4, hx4,sx4,kx4,pvalues4] = ACD_logburr(inx,[-0.12;0.13;0.52;1.35;1.35],2,[]);
    else
        [parameters4, stderrors4, LLF4, mux4, hx4,sx4,kx4,pvalues4] = ACD_logburr(inx,[],2,[]);% -316
    end 
    final_llf(item,20)=LLF4;
    params{item,20} = [parameters4,stderrors4,pvalues4];
    %%
    inz = insamplez(inx,mux4','bacd',parameters4);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux4(end),'logbacd2',parameters4);
    inzs{item}(:,20) = inz;
    outzs{item}(:,20) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(20,:)=[inp,in2p,outp,out2p,inpp,outpp];
    %%
    nu_len = 1;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    x0 = [alpha;beta;gamma;phi01];  
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(1,1)];
    upperbounds  = [1;1;1;10e8*ones(1,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood1(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood1(xx,inx,2);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood1(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux5]=logconstant_loglikelihood1(xx,inx,2);
        
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters21, stderrors21, LLF21, mux21, hx21, pvalues21,AIC21,BIC21] = ACD_logconstant_nnsnp1(inx,xx,2,[]) ;
    params{item,21} = [parameters21,stderrors21,pvalues21];
    final_llf(item,21)=llf;
    inz = insamplez(inx,mux5','led1',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux5(end),'logled12',xx);
    inzs{item}(:,21) = inz;
    outzs{item}(:,21) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(21,:)=[inp,in2p,outp,out2p,inpp,outpp];
    
    nu_len = 2;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02]; 
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(2,1)];
    upperbounds  = [1;1;1;10e8*ones(2,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood2(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood2(xx,inx,2);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood2(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux6]=logconstant_loglikelihood2(xx,inx,2);
        xx
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters22, stderrors22, LLF22, mux22, hx22,sx22,kx22, pvalues22,AIC22,BIC22] = ACD_logconstant_nnsnp(inx,xx,2,[]) ;
    params{item,22} = [parameters22,stderrors22,pvalues22];
    final_llf(item,22)=llf;
    inz = insamplez(inx,mux6','led2',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux6(end),'logled22',xx);
    inzs{item}(:,22) = inz;
    outzs{item}(:,22) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(22,:)=[inp,in2p,outp,out2p,inpp,outpp];
    nu_len = 3;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02;phi03]; 
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(3,1)];
    upperbounds  = [1;1;1;10e8*ones(3,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood3(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood3(xx,inx,2);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood3(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux7]=logconstant_loglikelihood3(xx,inx,2);
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters23, stderrors23, LLF23, mux23, hx23,sx23,kx23, pvalues23,AIC23,BIC23] = ACD_logconstant_nnsnp3(inx,xx,2,[]) ;
    params{item,23} = [parameters23,stderrors23,pvalues23];
    final_llf(item,23)=llf;
    inz = insamplez(inx,mux7','led3',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux7(end),'logled32',xx);
    inzs{item}(:,23) = inz;
    outzs{item}(:,23) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    %close outh
    [outc,outpp]=goodnesstest(outz,counts);
    test_results{item}(23,:)=[inp,in2p,outp,out2p,inpp,outpp];
    nu_len = 4;
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    phi04 = 0.05;
    x0 = [alpha;beta;gamma;phi01;phi02;phi03;phi04];    
    [r c] = size(x0);
    A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
    b = [zeros(1,3), 1 - 1e-8];
    lowerbounds  = [1e-8*ones(3,1);-10e8*ones(4,1)];
    upperbounds  = [1;1;1;10e8*ones(4,1)];
    rng(1)
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective',@(params)logconstant_loglikelihood4(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
    xx = run(gs,problem);
    llf=logconstant_loglikelihood4(xx,inx,2);
    
    %final_llf(1,nu_len)=llf;
    
    for i=2:20
        problem = createOptimProblem('fmincon','x0',xx, ...
            'objective',@(params)logconstant_loglikelihood4(params,inx,2),'Aineq',A,'bineq',b,'lb',lowerbounds,'ub',upperbounds);
        xx = run(gs,problem);
        [llf,mux8]=logconstant_loglikelihood4(xx,inx,2);
        parameters{nu_len}(i,:)=xx;
        %parameters{nu_len,i}=xx;
    end
    [parameters24, stderrors24, LLF24, mux24, hx24,sx24,kx24, pvalues24,AIC24,BIC24] = ACD_logconstant_nnsnp4(inx,xx,2,[]) ;
    params{item,24} = [parameters24,stderrors24,pvalues24];
    final_llf(item,24)=llf;
    inz = insamplez(inx,mux8','led4',xx);
    [inh,inp,~,~] = lbqtest(inz,Lags=15);
    [in2h,in2p,~,~] = lbqtest(inz.^2,Lags=15);
    
    outz = outsamplez(outx,outx1,mux8(end),'logled42',xx);
    inzs{item}(:,24) = inz;
    outzs{item}(:,24) = outz;
    [outh,outp,~,~] = lbqtest(outz,Lags=15);
    [out2h,out2p,~,~] = lbqtest(outz.^2,Lags=15);
    
    counts=histcounts(inz,20);
    %close inh1f
    [inc,inpp]=goodnesstest(inz,counts);
    
    outh=histogram(outz,20);
    counts=histcounts(outz,20);
    test_results{item}(24,:)=[inp,in2p,outp,out2p,inpp,outpp];



end
save('volume_all_acd_results.mat')

function [llf,mux] = constant_loglikelihood1(theta,x)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
psi = zeros(1,n);
psi(1) = alpha/(1-beta-gamma);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu_i_square = nu0^2+nu1^2;
a_nu_i = 1/(1-(2*nu1-2*nu1^2)/nu_i_square);
mu_z_2 = 2*(2*nu1^2)/nu_i_square-4*(2*nu1-2*nu1^2)/nu_i_square+2;
for i = 2:n
    psi(i) = alpha + beta*x(i-1) + gamma*psi(i-1);
    
    part0(i) = log((1+nu1*(-x(i)/psi(i)/a_nu_i+1))^2);
end
logfunction = - log(psi)-log(a_nu_i) - x'./(psi*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
llf = -ll;

mux = psi;
hx = psi.^2.*(a_nu_i.^2.*mu_z_2-1);
end
function [llf,mux] = constant_loglikelihood2(theta,x)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
phi02 = theta(5);
psi = zeros(1,n);
psi(1) = alpha/(1-beta-gamma);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu2 = phi02;
nu_i_square = nu0^2+nu1^2+nu2^2;
a_nu_i = 1/(1-(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square);
mu_z_2 = 2*(2*nu2+2*nu1^2-8*nu1*nu2+10*nu2^2)/nu_i_square...
    -4*(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square+2;
mu_z_3 = -6*(6*nu1*nu2-12*nu2^2)/nu_i_square+18*(2*nu2+2*nu1^2-8*nu1*nu2+10*nu2^2)/nu_i_square...
    -18*(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square+6;
mu_z_4 = 24*6*nu2^2/nu_i_square...
    -96*(6*nu1*nu2-12*nu2^2)/nu_i_square...
    +144*(2*nu2+2*nu1^2-8*nu1*nu2+10*nu2^2)/nu_i_square...
    -96*(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square+24;
for i = 2:n
    psi(i) = alpha + beta*x(i-1) + gamma*psi(i-1);
    
    part0(i) = log((1+nu1*(-x(i)/psi(i)/a_nu_i+1)+nu2*L_2(x(i)/psi(i)/a_nu_i))^2);
end
logfunction = - log(psi)-log(a_nu_i) - x'./(psi*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
llf = -ll;

mux = psi;
hx = psi.^2.*(a_nu_i.^2.*mu_z_2-1);
end


function [llf,mux] = constant_loglikelihood3(theta,x)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
phi02 = theta(5);
phi03 = theta(6);
psi = zeros(1,n);
psi(1) = alpha/(1-beta-gamma);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu2 = phi02;
nu3 = phi03;
nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2;
a_nu_i = 1/(1-(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square);
mu_z_2 = 2*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+2*nu1^2+10*nu2^2+24*nu3^2)/nu_i_square...
    -4*(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square+2;
mu_z_3 = -6*(2*nu3+6*nu1*nu2-12*nu1*nu3+48*nu2*nu3-12*nu2^2-56*nu3^2)/nu_i_square...
    +18*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+2*nu1^2+10*nu2^2+24*nu3^2)/nu_i_square...
    -18*(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square...
    +6;
mu_z_4 = 24*(8*nu1*nu3-48*nu2*nu3+6*nu2^2+78*nu3^2)/nu_i_square...
    -96*(2*nu3+6*nu1*nu2-12*nu1*nu3+48*nu2*nu3-12*nu2^2-56*nu3^2)/nu_i_square...
    +144*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+2*nu1^2+10*nu2^2+24*nu3^2)/nu_i_square...
    -96*(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square...
    +24;
for i = 2:n
    psi(i) = alpha + beta*x(i-1) + gamma*psi(i-1);
    
    part0(i) = log((1+nu1*(-x(i)/psi(i)/a_nu_i+1)+nu2*L_2(x(i)/psi(i)/a_nu_i)+nu3*L_3(x(i)/psi(i)/a_nu_i))^2);
end
logfunction = - log(psi)-log(a_nu_i) - x'./(psi*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
llf = -ll;

mux = psi;
hx = psi.^2.*(a_nu_i.^2.*mu_z_2-1);
sx = (a_nu_i.^3.*mu_z_3+2-3*a_nu_i.^2.*mu_z_2)./((a_nu_i.^2.*mu_z_2-1)).^(3/2);
kx = (a_nu_i.^4.*mu_z_4-3+4*a_nu_i.^2.*mu_z_2-4*a_nu_i.^3.*mu_z_3)./((a_nu_i.^2.*mu_z_2-1)).^2;
end

function [llf,mux] = constant_loglikelihood4(theta,x)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
phi02 = theta(5);
phi03 = theta(6);
phi04 = theta(7);
psi = zeros(1,n);
psi(1) = alpha/(1-beta-gamma);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu2 = phi02;
nu3 = phi03;
nu4 = phi04;
nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2+nu4^2;
a_nu_i = 1/(1-(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square);
mu_z_2 = 2*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square...
    -4*(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square+2;
mu_z_3 = -6*(2*nu3+6*nu1*nu2-12*nu1*nu3+8*nu1*nu4+48*nu2*nu3-48*nu2*nu4+156*nu3*nu4-12*nu2^2-56*nu3^2-152*nu4^2)/nu_i_square...
    +18*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square...
    -18*(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square...
    +6;
mu_z_4 = 24*(2*nu4+8*nu1*nu3-16*nu1*nu4-48*nu2*nu3+88*nu2*nu4-304*nu3*nu4+6*nu2^2+78*nu3^2+346*nu4^2)/nu_i_square...
    -96*(2*nu3+6*nu1*nu2-12*nu1*nu3+8*nu1*nu4+48*nu2*nu3-48*nu2*nu4+156*nu3*nu4-12*nu2^2-56*nu3^2-152*nu4^2)/nu_i_square...
    +144*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square...
    -96*(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square...
    +24;
for i = 2:n
    psi(i) = alpha + beta*x(i-1) + gamma*psi(i-1);
    
    part0(i) = log((1+nu1*(-x(i)/psi(i)/a_nu_i+1)+nu2*L_2(x(i)/psi(i)/a_nu_i)+nu3*L_3(x(i)/psi(i)/a_nu_i)+nu4*L_4(x(i)/psi(i)/a_nu_i))^2);
end
logfunction = - log(psi)-log(a_nu_i) - x'./(psi*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
llf = -ll;

mux = psi;
hx = psi.^2.*(a_nu_i.^2.*mu_z_2-1);
sx = (a_nu_i.^3.*mu_z_3+2-3*a_nu_i.^2.*mu_z_2)./((a_nu_i.^2.*mu_z_2-1)).^(3/2);
kx = (a_nu_i.^4.*mu_z_4-3+4*a_nu_i.^2.*mu_z_2-4*a_nu_i.^3.*mu_z_3)./((a_nu_i.^2.*mu_z_2-1)).^2;
end
function [ll,mux] = logconstant_loglikelihood1(theta,x,type)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
psi = zeros(1,n);
psi(1) = mean(x);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu_i_square = nu0^2+nu1^2;
a_nu_i = 1/(1-(2*nu1-2*nu1^2)/nu_i_square);
mu_z_2 = 2*(2*nu1^2)/nu_i_square-4*(2*nu1-2*nu1^2)/nu_i_square+2;
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/exp(psi(1));
for i = 2:n
    if type==1
        psi(i) = alpha + beta*log(x(i-1)) + gamma*psi(i-1);
    else
        
        psi(i) = alpha + beta*epsilon(i-1) + gamma*psi(i-1);
        epsilon(i) = x(i)/exp(psi(i));
    end
    
    part0(i) = log((1+nu1*(-x(i)/exp(psi(i))/a_nu_i+1))^2);
end
logfunction = - psi-log(a_nu_i) - x'./(exp(psi)*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
ll = -ll;

mux = exp(psi);
hx = exp(psi).^2.*(a_nu_i.^2.*mu_z_2-1);
end
function [ll,mux] = logconstant_loglikelihood2(theta,x,type)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
phi02 = theta(5);
psi = zeros(1,n);
psi(1) = mean(x);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu2 = phi02;
nu_i_square = nu0^2+nu1^2+nu2^2;
a_nu_i = 1/(1-(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square);
mu_z_2 = 2*(2*nu2+2*nu1^2-8*nu1*nu2+10*nu2^2)/nu_i_square...
    -4*(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square+2;
mu_z_3 = -6*(6*nu1*nu2-12*nu2^2)/nu_i_square+18*(2*nu2+2*nu1^2-8*nu1*nu2+10*nu2^2)/nu_i_square...
    -18*(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square+6;
mu_z_4 = 24*6*nu2^2/nu_i_square...
    -96*(6*nu1*nu2-12*nu2^2)/nu_i_square...
    +144*(2*nu2+2*nu1^2-8*nu1*nu2+10*nu2^2)/nu_i_square...
    -96*(2*nu1-2*nu1^2+4*nu1*nu2-4*nu2^2)/nu_i_square+24;
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/exp(psi(1));
for i = 2:n
    if type==1
        psi(i) = alpha + beta*log(x(i-1)) + gamma*psi(i-1);
    else
        
        psi(i) = alpha + beta*epsilon(i-1) + gamma*psi(i-1);
        epsilon(i) = x(i)/exp(psi(i));
    end
    
    part0(i) = log((1+nu1*(-x(i)/exp(psi(i))/a_nu_i+1)+nu2*L_2(x(i)/exp(psi(i))/a_nu_i))^2);
end
logfunction = - psi-log(a_nu_i) - x'./(exp(psi)*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
ll = -ll;

mux = exp(psi);
hx = exp(psi).^2.*(a_nu_i.^2.*mu_z_2-1);
sx = (a_nu_i.^3.*mu_z_3+2-3*a_nu_i.^2.*mu_z_2)./((a_nu_i.^2.*mu_z_2-1)).^(3/2);
kx = (a_nu_i.^4.*mu_z_4-3+4*a_nu_i.^2.*mu_z_2-4*a_nu_i.^3.*mu_z_3)./((a_nu_i.^2.*mu_z_2-1)).^2;
end


function [ll,mux] = logconstant_loglikelihood3(theta,x,type)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
phi02 = theta(5);
phi03 = theta(6);
psi = zeros(1,n);
psi(1) = mean(x);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu2 = phi02;
nu3 = phi03;
nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2;
a_nu_i = 1/(1-(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square);
mu_z_2 = 2*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+2*nu1^2+10*nu2^2+24*nu3^2)/nu_i_square...
    -4*(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square+2;
mu_z_3 = -6*(2*nu3+6*nu1*nu2-12*nu1*nu3+48*nu2*nu3-12*nu2^2-56*nu3^2)/nu_i_square...
    +18*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+2*nu1^2+10*nu2^2+24*nu3^2)/nu_i_square...
    -18*(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square...
    +6;
mu_z_4 = 24*(8*nu1*nu3-48*nu2*nu3+6*nu2^2+78*nu3^2)/nu_i_square...
    -96*(2*nu3+6*nu1*nu2-12*nu1*nu3+48*nu2*nu3-12*nu2^2-56*nu3^2)/nu_i_square...
    +144*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+2*nu1^2+10*nu2^2+24*nu3^2)/nu_i_square...
    -96*(2*nu1+4*nu1*nu2+6*nu2*nu3-2*nu1^2-4*nu2^2-6*nu3^2)/nu_i_square...
    +24;
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/exp(psi(1));
for i = 2:n
    if type==1
        psi(i) = alpha + beta*log(x(i-1)) + gamma*psi(i-1);
    else
        
        psi(i) = alpha + beta*epsilon(i-1) + gamma*psi(i-1);
        epsilon(i) = x(i)/exp(psi(i));
    end
    
    part0(i) = log((1+nu1*(-x(i)/exp(psi(i))/a_nu_i+1)+nu2*L_2(x(i)/exp(psi(i))/a_nu_i)+nu3*L_3(x(i)/exp(psi(i))/a_nu_i))^2);
end
logfunction = - psi-log(a_nu_i) - x'./(exp(psi)*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
ll = -ll;

mux = exp(psi);
hx = exp(psi).^2.*(a_nu_i.^2.*mu_z_2-1);
sx = (a_nu_i.^3.*mu_z_3+2-3*a_nu_i.^2.*mu_z_2)./((a_nu_i.^2.*mu_z_2-1)).^(3/2);
kx = (a_nu_i.^4.*mu_z_4-3+4*a_nu_i.^2.*mu_z_2-4*a_nu_i.^3.*mu_z_3)./((a_nu_i.^2.*mu_z_2-1)).^2;
end


function [ll,mux] = logconstant_loglikelihood4(theta,x,type)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta = theta(2);
gamma = theta(3);
phi01 = theta(4);
phi02 = theta(5);
phi03 = theta(6);
phi04 = theta(7);
psi = zeros(1,n);
psi(1) = mean(x);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu2 = phi02;
nu3 = phi03;
nu4 = phi04;
nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2+nu4^2;
a_nu_i = 1/(1-(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square);
mu_z_2 = 2*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square...
    -4*(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square+2;
mu_z_3 = -6*(2*nu3+6*nu1*nu2-12*nu1*nu3+8*nu1*nu4+48*nu2*nu3-48*nu2*nu4+156*nu3*nu4-12*nu2^2-56*nu3^2-152*nu4^2)/nu_i_square...
    +18*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square...
    -18*(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square...
    +6;
mu_z_4 = 24*(2*nu4+8*nu1*nu3-16*nu1*nu4-48*nu2*nu3+88*nu2*nu4-304*nu3*nu4+6*nu2^2+78*nu3^2+346*nu4^2)/nu_i_square...
    -96*(2*nu3+6*nu1*nu2-12*nu1*nu3+8*nu1*nu4+48*nu2*nu3-48*nu2*nu4+156*nu3*nu4-12*nu2^2-56*nu3^2-152*nu4^2)/nu_i_square...
    +144*(2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square...
    -96*(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square...
    +24;
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/exp(psi(1));
for i = 2:n
    if type==1
        psi(i) = alpha + beta*log(x(i-1)) + gamma*psi(i-1);
    else
        
        psi(i) = alpha + beta*epsilon(i-1) + gamma*psi(i-1);
        epsilon(i) = x(i)/exp(psi(i));
    end
    
    part0(i) = log((1+nu1*(-x(i)/exp(psi(i))/a_nu_i+1)+nu2*L_2(x(i)/exp(psi(i))/a_nu_i)+nu3*L_3(x(i)/exp(psi(i))/a_nu_i)+nu4*L_4(x(i)/exp(psi(i))/a_nu_i))^2);
end
logfunction = - psi-log(a_nu_i) - x'./(exp(psi)*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
ll = -ll;

mux = exp(psi);
hx = exp(psi).^2.*(a_nu_i.^2.*mu_z_2-1);
sx = (a_nu_i.^3.*mu_z_3+2-3*a_nu_i.^2.*mu_z_2)./((a_nu_i.^2.*mu_z_2-1)).^(3/2);
kx = (a_nu_i.^4.*mu_z_4-3+4*a_nu_i.^2.*mu_z_2-4*a_nu_i.^3.*mu_z_3)./((a_nu_i.^2.*mu_z_2-1)).^2;
end


function re = L_2(x)
re = 1/2*(x^2-4*x+2);
end
function re = L_3(x)
re = 1/6*(-x^3+9*x^2-18*x+6);
end
function re = L_4(x)
re = 1/24*(x^4-16*x^3+72*x^2-96*x+24);
end