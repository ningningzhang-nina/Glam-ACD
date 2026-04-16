function [ll,loglikelihood,mux,hx,sx,kx] = loglikelihood_logburr(parameters,x,type)
n = length(x);
alpha = parameters(1);
beta = parameters(2);
Gamma = parameters(3);
kappa = parameters(4);
sigma2 = parameters(5);
part = gamma(1+1/kappa)*gamma(1/sigma2-1/kappa)/((sigma2)^(1+1/kappa)*gamma(1/sigma2+1));
xi = zeros(1,n);
xi(1) = alpha/(1-beta-Gamma)/part;
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/(exp(xi(1))*part);
for i=2:length(x)
    if type==1
        xi(i) = alpha/part + beta/part*log(x(i-1)) + Gamma*xi(i-1);
    else
        
        xi(i) = alpha/part + beta/part*epsilon(i-1) + Gamma*xi(i-1);
        epsilon(i) = x(i)/(exp(xi(i))*part);
    end
end
loglikelihood = log(kappa)-kappa*xi+(kappa-1)*log(x')...
    -(1/sigma2+1)*log(1+sigma2*exp(xi).^(-kappa).*x'.^(kappa));
ll = sum(loglikelihood);
ll = -ll;
mux = exp(xi).*part;
hx = mux.^2;
sx = 2;
kx = 6;
end