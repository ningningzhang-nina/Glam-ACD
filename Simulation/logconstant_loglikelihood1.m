function [ll,logfunction,mux,hx] = logconstant_loglikelihood1(theta,x,type)

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

