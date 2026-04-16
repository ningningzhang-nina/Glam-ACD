function [ll,loglikelihood,mux,hx,sx,kx] = loglikelihood_generalized_gamma(parameters,x)
Alpha = parameters(1);
beta = parameters(2);
Gamma = parameters(3);
delta = parameters(4);
alpha = parameters(5);
theta = zeros(1,length(x));
kappa_alpha_delta = gamma(alpha+1/delta)/gamma(alpha);
theta(1) = Alpha/(1-beta-Gamma)/kappa_alpha_delta;
for i=2:length(x)
    theta(i) = Alpha/kappa_alpha_delta + beta/kappa_alpha_delta*x(i-1) + Gamma*theta(i-1);
    %theta(i) = kappa_alpha_delta^(Gamma-1)*exp(Gamma*log(theta(i-1))+beta/kappa_alpha_delta*x(i-1)/theta(i-1)+Alpha);
end
loglikelihood = log(delta)+(delta*alpha-1)*log(x')-delta*alpha*log(theta)-log(gamma(alpha))-(x'./theta).^delta;
ll = sum(loglikelihood);
ll = -ll;
mux = theta*kappa_alpha_delta;
hx = theta.^2*(gamma(alpha+2/delta)/gamma(alpha)-(gamma(alpha+1/delta)/gamma(alpha))^2);
ex2 = theta.^2*gamma(alpha+2/delta)/gamma(alpha);
ex3 = theta.^3*gamma(alpha+3/delta)/gamma(alpha);
ex4 = theta.^4*gamma(alpha+4/delta)/gamma(alpha);
sx = (ex3+2*mux.^3-3*ex2.*mux)./(hx).^(3/2);
kx = (ex4-3*mux.^4+4*ex2.*mux.^2-4*ex3.*mux)./(hx).^(2);
% pdTrue = GeneralizedGamma(1/theta(2), delta, alpha);
% x = 1:1:500;
% f = pdTrue.pdf(x);
% plot(x, f);
end