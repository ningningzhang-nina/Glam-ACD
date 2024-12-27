function [f,mu,h,s,k]=compute_moments(x,model,parameters)
    if model == "eacd"
        lambda = parameters;
        f=lambda.*exp(-lambda.*x);
        mu=1/lambda;
        h=1/(lambda^2);
        s = 2;
        k = 6;
    elseif model =="wacd"
        lambda=parameters(1);
        k = parameters(2);
        f=k/lambda*(x./lambda).^(k-1).*exp(-(x./lambda).^k);
        mu1 = mur_weibull(lambda,k,1);
        mu2 = mur_weibull(lambda,k,2);
        mu3 = mur_weibull(lambda,k,3);
        mu4 = mur_weibull(lambda,k,4);
        [h,s,k]=compute_high_moments(mu1,mu2,mu3,mu4);
        mu=mu1;
    elseif model=="ggacd"
        a=parameters(1);
        delta = parameters(2);
        kappa = parameters(3);
        f=(delta/(a^(delta*kappa))).*x.^(delta*kappa-1).*exp(-(x./a).^delta)./gamma(kappa);
        mu1 = mur_gg(a,delta,kappa,1);
        mu2 = mur_gg(a,delta,kappa,2);
        mu3 = mur_gg(a,delta,kappa,3);
        mu4 = mur_gg(a,delta,kappa,4);
        mu=mu1;
        [h,s,k]=compute_high_moments(mu1,mu2,mu3,mu4);
    elseif model=="bacd"
        lambda=parameters(1);
        kappa = parameters(2);
        sigma2 = parameters(3);
        f=kappa/sigma2/lambda.*(x./lambda).^(kappa-1).*(1+(x./lambda).^(kappa)).^(-1/sigma2-1);
        mu1 = mur_burr(lambda,sigma2,kappa,1);
        mu2 = mur_burr(lambda,sigma2,kappa,2);
        mu3 = mur_burr(lambda,sigma2,kappa,3);
        mu4 = mur_burr(lambda,sigma2,kappa,4);
        mu=mu1;
        [h,s,k]=compute_high_moments(mu1,mu2,mu3,mu4);
    elseif model=="led1"
        nu0 = 1;
        nu1 = parameters(1);
        [mu,mu2,mu3,mu4] = mur_led(nu0,nu1,0,0,0);
        [h,s,k]=compute_high_moments(mu,mu2,mu3,mu4);
        f=density_led(x,nu1,0,0,0);
    elseif model=="led2"
        nu0 = 1;
        nu1 = parameters(1);
        nu2 = parameters(2);
        [mu,mu2,mu3,mu4] = mur_led(nu0,nu1,nu2,0,0);
        [h,s,k]=compute_high_moments(mu,mu2,mu3,mu4);
        f=density_led(x,nu1,nu2,0,0);
    elseif model=="led3"
        nu0 = 1;
        nu1 = parameters(1);
        nu2 = parameters(2);
        nu3 = parameters(3);
        [mu,mu2,mu3,mu4] = mur_led(nu0,nu1,nu2,nu3,0);
        [h,s,k]=compute_high_moments(mu,mu2,mu3,mu4);
        f=density_led(x,nu1,nu2,nu3,0);
    else
        nu0 = 1;
        nu1 = parameters(1);
        nu2 = parameters(2);
        nu3 = parameters(3);
        nu4 = parameters(4);
        [mu,mu2,mu3,mu4] = mur_led(nu0,nu1,nu2,nu3,nu4);
        [h,s,k]=compute_high_moments(mu,mu2,mu3,mu4);
        f=density_led(x,nu1,nu2,nu3,nu4);
    end
end
function mur = mur_weibull(lambda,k,r)
mur = lambda.^r*gamma(1+r/k);
end
function mur = mur_gg(a,delta,kappa,r)
mur = a.^r*gamma(kappa+r/delta)/gamma(kappa);
end
function mur = mur_burr(lambda,sigma2,kappa,r)
mur = lambda.^r*gamma(1+r/kappa)*gamma(1/sigma2-r/kappa)/(sigma2)/gamma(1/sigma2+1);
end
function [mu1,mu2,mu3,mu4] = mur_led(nu0,nu1,nu2,nu3,nu4)
nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2+nu4^2;
gamma_1 = (2*nu0*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu4^2)/nu_i_square;
gamma_2 = (2*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+10*nu2^2+24*nu3^2+44*nu4^2)/nu_i_square;
gamma_3 = (2*nu3+6*nu1*nu2-12*nu1*nu3+8*nu1*nu4+48*nu2*nu3-48*nu2*nu4+156*nu3*nu4-12*nu2^2-56*nu3^2-152*nu4^2)/nu_i_square;
gamma_4 = (2*nu4+8*nu1*nu3-16*nu1*nu4-48*nu2*nu3+88*nu2*nu4-304*nu3*nu4+6*nu2^2+78*nu3^2+346*nu4^2)/nu_i_square;

mu1 = 1-gamma_1;

%a_nu_i = 1/(1-(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square);
mu2 = 2*gamma_2-4*gamma_1+2;
mu3 = -6*gamma_3+18*gamma_2-18*gamma_1+6;
mu4 = 24*gamma_4-96*gamma_3+144*gamma_2-96*gamma_1+24;
end
function f=density_led(x,nu1,nu2,nu3,nu4)
    f=exp(-x)./(1+nu1.^2+nu2.^2+nu3.^2+nu4.^2).*(Lague(x,0)+nu1*Lague(x,1)+nu2*Lague(x,2)+nu3*Lague(x,3)+nu4*Lague(x,4)).^2;
end
function l=Lague(x,k)
   if k==0
       l=1;
   elseif k==1
       l=-x+1;
   elseif k==2
       l=(x.^2-4.*x+2)./2;
   elseif k==3
       l=(-x.^3+9.*x.^2-18.*x+6)./6;
   else
       l=(x.^4-16.*x.^3+72.*x.^2-96.*x+24)./24;
   end
end
function [hx,sx,kx]=compute_high_moments(mu1,mu2,mu3,mu4)
hx = -mu1.^2+mu2;
sx=(2*mu1.^3-3*mu1*mu2+mu3)/(hx^(3/2));
kx=(-3*mu1.^4+6*mu1.^2.*mu2-4*mu1.*mu3+mu4)/(hx^(2));
end