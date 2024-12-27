function [parameters, stderrors, LLF, mux, hx, sx, kx, pvalues,AIC,BIC] = ACD_constant_nnsnp3(x,startingvalues,options) 
% Specify Initial Values
if isempty(startingvalues)
    alpha = 0.0084;
    beta = 0.0405;
    gamma = 0.9511;
    phi01 = -0.67;
    phi02 = 0.05;
    phi03 = 0.05;
    startingvalues = [alpha;beta;gamma;phi01;phi02;phi03];    
 end 

% Specify constraints used by fmincon (A, b) and lower and upper bounds of parameters
% Example: alpha>0,beta>0,gamma>0,beta+gamma<1,
[r c] = size(startingvalues);
A = [-eye(3) zeros(3,r-3);0 ones(1,2) zeros(1,r-3)]; 
%; 
A
b = [zeros(1,3), 1 - 1e-8];
% 
b
% lowerbounds  = [1e-8*ones(3,1);-1;1e-8;-1*ones(4,1);1e-8;-1*ones(3,1)];
% lowerbounds
% upperbounds  = [1;1;1;1;1;ones(4,1);1;ones(3,1)];
% upperbounds
lowerbounds  = [1e-8*ones(3,1);-10e8*ones(3,1)];
lowerbounds
upperbounds  = [1;1;1;10e8*ones(3,1)];
upperbounds
% define options
if isempty(options)
    options  =  optimset('fmincon');
    options  =  optimset(options , 'Algorithm ','sqp');
    options  =  optimset(options , 'TolFun'      , 1e-006);
    options  =  optimset(options , 'TolX'        , 1e-006);
    options  =  optimset(options , 'TolCon'      , 1e-006);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxIter'     , 1500);
    options  =  optimset(options , 'Jacobian'     ,'off');
    options  =  optimset(options , 'MeritFunction'     ,'multiobj');
    options  =  optimset(options , 'MaxFunEvals' , 6000);
end
%'myconfun'
[parameters, LLF, EXITFLAG, OUTPUT, Lambda, GRAD, HESSIAN] =  fmincon('constant_loglikelihood3', startingvalues ,A, b, [],[],lowerbounds, upperbounds, [], options, x);

if EXITFLAG<=0
   EXITFLAG
   fprintf(1,'Convergence has been not successful!\n')
end

[LLF,likelihoods,mux,hx,sx,kx] = constant_loglikelihood3(parameters, x);
LLF = -LLF;


HESSIAN = zeros(r,r);
epsilon = 0.00000003 * parameters;
for i=1:r
    for j = 1:r
      x1 = parameters;
      x2 = parameters;
      x3 = parameters;
      x4 = parameters;
      x1(i) = x1(i) + epsilon(i); 
      x1(j) = x1(j) + epsilon(j);
      x2(i) = x2(i) + epsilon(i); 
      x2(j) = x2(j) - epsilon(j);
      x3(i) = x3(i) - epsilon(i); 
      x3(j) = x3(j) + epsilon(j);
      x4(i) = x4(i) - epsilon(i); 
      x4(j) = x4(j) - epsilon(j);
      [LLF1,~,~, ~, ~,~] = constant_loglikelihood3(x1, x);
      [LLF2,~,~, ~, ~,~] = constant_loglikelihood3(x2, x);
      [LLF3,~,~, ~, ~,~] = constant_loglikelihood3(x3, x);
      [LLF4,~,~, ~, ~,~] = constant_loglikelihood3(x4, x);
      HESSIAN(i, j) = (LLF1-LLF2-LLF3+LLF4)/(4*epsilon(i)*epsilon(j));
    end
end
% 
% HESSIAN
% asymptotic standard errors
stderrors = sqrt(diag(inv(HESSIAN))); 
stderrors
% t-statistics
tstats = parameters./stderrors;
tstats
% robust standard errors
h=parameters*eps;
hplus=parameters+h;
hminus=parameters-h;
T = length(x);
m = 0;
likelihoodsplus=zeros(T-m,length(parameters));
likelihoodsminus=zeros(T-m,length(parameters));
for i=1:length(parameters)
    hparameters=parameters;
    hparameters(i)=hplus(i);
    [holder1,indivlike] =  constant_loglikelihood3(hparameters, x);
    likelihoodsplus(:,i)=indivlike';
end
for i=1:length(parameters)
    hparameters=parameters;
    hparameters(i)=hminus(i);
     [holder1, indivlike] =  constant_loglikelihood3(hparameters, x);
    likelihoodsminus(:,i)=indivlike';
end
scores=(likelihoodsplus-likelihoodsminus)./(2*repmat(h',T-m,1));
scores=scores-repmat(mean(scores),T-m,1);
S=scores'*scores;
robustSE = diag((HESSIAN^(-1))*S*(HESSIAN^(-1)));
robusttstats = parameters./robustSE;
pvalues = 2*(1-normcdf(abs(tstats)));

% Saving and Organizing Results
summary.Iterations = OUTPUT.iterations;

% Statistics to be saved
A = char('Coeff', 'SErrors', 'Tstats', 'pvalues','RobustSErrors', 'RobustTstats');

% Result vectors
C = char('parameters', 'stderrors', 'tstats','pvalues', 'robustSE', 'robusttstats');

% Parameters to be saved   
B = strcat(char('alpha','beta','gamma','phi01','phi02','phi03'));
   
for i = 1:size(A,1);
    for j = 1:size(B,1);
        eval(['summary.',strcat(A(i,:)),'.',strcat(B(j,:)),' =',strcat(C(i,:)),'(j);']);
    end
    clear j
end
clear i

summary.Scores=scores;
summary.LLF = LLF;
AIC = -2*LLF+2*size(parameters,1);
BIC = -2*LLF+size(parameters,1)*log(length(x));
summary.AIC = -2*LLF+2*size(parameters,1);
summary.BIC = -2*LLF+size(parameters,1)*log(length(x));
fprintf('------------------------------------------------------------\n')
fprintf('Convergence achieved after %1.0f iterations\n', OUTPUT.iterations)
fprintf('------------------------------------------------------------\n')
fprintf('Parameters  Coefficients  Std Errors    T-stats      pvalues\n')
fprintf('------------------------------------------------------------\n')
for i = 1:size(parameters,1)
    if parameters(i) < 0
        fprintf(strcat('  %s     %1.',num2str(round(5-length(sprintf('%1.0f', abs(parameters(i)))))),'f      %1.',num2str(round(5-length(sprintf('%1.0f', stderrors(i))))),'f       %1.',num2str(round(5-length(sprintf('%1.0f',abs(tstats(i)))))),'f       %1.',num2str(round(5-length(sprintf('%1.0f', pvalues(i))))),'f\n'), B(i,:), parameters(i), stderrors(i), tstats(i),pvalues(i))
    else    
        fprintf(strcat('  %s      %1.',num2str(round(5-length(sprintf('%1.0f', abs(parameters(i)))))),'f      %1.',num2str(round(5-length(sprintf('%1.0f', stderrors(i))))),'f       %1.',num2str(round(5-length(sprintf('%1.0f', abs(tstats(i)))))),'f       %1.',num2str(round(5-length(sprintf('%1.0f', pvalues(i))))),'f\n'), B(i,:), parameters(i), stderrors(i), tstats(i),pvalues(i))
    end
end
fprintf('------------------------------------------------------------\n')
fprintf('Log Likelihood: %1.0f\n', LLF)
fprintf('Akaike Information Criteron: %1.0f\n', summary.AIC)
fprintf('Bayesian Information Criteron: %1.0f\n', summary.BIC)
fprintf('------------------------------------------------------------\n')

end

