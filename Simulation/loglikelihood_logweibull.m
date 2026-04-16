function [ll,loglikelihood,mux,hx,sx,kx] = loglikelihood_logweibull(parameters,x,type)
alpha = parameters(1);
beta = parameters(2);
Gamma = parameters(3);
k = parameters(4);
psi = zeros(1,length(x));
psi(1) = alpha/(1-beta-Gamma);
epsilon = zeros(1,length(x));
epsilon(1) = x(1)/exp(psi(1));
for i=2:length(x)
    if type==1
        psi(i) = alpha + beta*log(x(i-1)) + Gamma*psi(i-1);
    else
        
        psi(i) = alpha + beta*epsilon(i-1) + Gamma*psi(i-1);
        epsilon(i) = x(i)/exp(psi(i));
    end
end
loglikelihood = log(k./x'.*((x'.*gamma(1+1/k))./exp(psi)).^k.*exp(-(((x'.*(gamma(1+1/k))./exp(psi)).^k))));
ll = sum(loglikelihood);
ll = -ll;
mux = exp(psi);
hx = exp(psi).^2/(gamma(1+1/k))^2*(gamma(1+2/k)-gamma(1+1/k))^2;
sx = (2*(gamma(1+1/k))^3-3*gamma(1+1/k)*gamma(1+2/k)+gamma(1+3/k))/(gamma(1+2/k)-(gamma(1+1/k))^2)^(3/2);
kx = (-6*(gamma(1+1/k))^4+12*(gamma(1+1/k))^2*gamma(1+2/k)-3*(gamma(1+2/k))^2-4*gamma(1+1/k)*gamma(1+3/k)+gamma(1+4/k))/(gamma(1+2/k)-(gamma(1+1/k))^2)^(2)+3;
end