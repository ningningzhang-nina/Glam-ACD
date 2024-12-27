function [ll,logfunction,mux,hx,sx,kx] = constant22_loglikelihood2(theta,x)

n = length(x);
nu0 = 1;
alpha = theta(1);
beta1 = theta(2);
gamma1 = theta(3);
beta2 = theta(4);
gamma2 = theta(5);
phi01 = theta(6);
phi02 = theta(7);
psi = zeros(1,n);
psi(1) = alpha/(1-beta1-gamma1-beta2-gamma2);
psi(2) = alpha + beta1*x(1) + gamma1*psi(1);
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
for i = 3:n
    psi(i) = alpha + beta1*x(i-1) + gamma1*psi(i-1)+ beta2*x(i-2) + gamma2*psi(i-2);
    
    part0(i) = log((1+nu1*(-x(i)/psi(i)/a_nu_i+1)+nu2*L_2(x(i)/psi(i)/a_nu_i))^2);
end
logfunction = - log(psi)-log(a_nu_i) - x'./(psi*a_nu_i) - log(nu_i_square)+part0;
ll = sum(logfunction);
ll = -ll;

mux = psi;
hx = psi.^2.*(a_nu_i.^2.*mu_z_2-1);
sx = (a_nu_i.^3.*mu_z_3+2-3*a_nu_i.^2.*mu_z_2)./((a_nu_i.^2.*mu_z_2-1)).^(3/2);
kx = (a_nu_i.^4.*mu_z_4-3+4*a_nu_i.^2.*mu_z_2-4*a_nu_i.^3.*mu_z_3)./((a_nu_i.^2.*mu_z_2-1)).^2;
end



function re = L_2(x)
re = 1/2*(x^2-4*x+2);
end