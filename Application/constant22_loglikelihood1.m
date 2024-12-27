function [ll,logfunction,mux,hx] = constant22_loglikelihood1(theta,x)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta1 = theta(2);
gamma1 = theta(3);
beta2 = theta(4);
gamma2 = theta(5);
phi01 = theta(6);
psi = zeros(1,n);
psi(1) = alpha/(1-beta1-gamma1-beta2-gamma2);
psi(2) = alpha + beta1*x(1) + gamma1*psi(1);
part0 = zeros(1,n);
part0(1) = 1;
nu1 = phi01;
nu_i_square = nu0^2+nu1^2;
a_nu_i = 1/(1-(2*nu1-2*nu1^2)/nu_i_square);
mu_z_2 = 2*(2*nu1^2)/nu_i_square-4*(2*nu1-2*nu1^2)/nu_i_square+2;
for i = 3:n
    psi(i) = alpha + beta1*x(i-1) + gamma1*psi(i-1) + beta2*x(i-2) + gamma2*psi(i-2);
    
    part0(i) = log((1+nu1*(-x(i)/psi(i)/a_nu_i+1))^2);
end
logfunction = - log(psi)-log(a_nu_i) - x'./(psi*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
ll = -ll;

mux = psi;
hx = psi.^2.*(a_nu_i.^2.*mu_z_2-1);
end

