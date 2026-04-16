function [ll,loglikelihood,mux,hx,sx,kx] = loglikelihood_burr(parameters,x)
n = length(x);
alpha = parameters(1);
beta = parameters(2);
Gamma = parameters(3);
kappa = parameters(4);
sigma2 = parameters(5);
part = gamma(1+1/kappa)*gamma(1/sigma2-1/kappa)/((sigma2)^(1+1/kappa)*gamma(1/sigma2+1));
xi = zeros(1,n);
xi(1) = alpha/(1-beta-Gamma)/part;
for i=2:n
    xi(i) = alpha/part + beta/part*x(i-1) + Gamma*xi(i-1);
end
loglikelihood = log(kappa)-kappa*log(xi)+(kappa-1)*log(x')...
    -(1/sigma2+1)*log(1+sigma2*xi.^(-kappa).*x'.^(kappa));
ll = sum(loglikelihood);
ll = -ll;
mux = xi.*part;
hx = mux.^2;
sx = 2;
kx = 6;
end