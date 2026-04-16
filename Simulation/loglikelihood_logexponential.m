function [ll,loglikelihood,mux,hx,sx,kx] = loglikelihood_logexponential(parameters,x,type)
alpha = parameters(1);
beta = parameters(2);
gamma = parameters(3);
psi = zeros(1,length(x));
psi(1) = alpha/(1-beta-gamma);
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/exp(psi(1));
for i=2:length(x)
    if type==1
        psi(i) = alpha + beta*log(x(i-1)) + gamma*psi(i-1);
    else
        
        psi(i) = alpha + beta*epsilon(i-1) + gamma*psi(i-1);
        epsilon(i) = x(i)/exp(psi(i));
    end
end
loglikelihood = -psi-x'./exp(psi);
ll = sum(loglikelihood);
ll = -ll;
mux = exp(psi);
hx = exp(psi).^2;
sx = 2;
kx = 6;
end