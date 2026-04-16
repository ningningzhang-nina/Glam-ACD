function [ll,loglikelihood,mux,hx,sx,kx] = loglikelihood_exponential(parameters,x)
alpha = parameters(1);
beta = parameters(2);
gamma = parameters(3);
psi = zeros(1,length(x));
psi(1) = alpha/(1-beta-gamma);
for i=2:length(x)
    psi(i) = alpha + beta*x(i-1) + gamma*psi(i-1);
end
loglikelihood = -log(psi)-x'./psi;
ll = sum(loglikelihood);
ll = -ll;
mux = psi;
hx = psi.^2;
sx = 2;
kx = 6;
end