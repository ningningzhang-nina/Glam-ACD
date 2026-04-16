function z = insamplez(xi,psi,model,parameters)
    
    if model=="eacd"
        
        z = 1-exp(-xi./psi);
    elseif model=="wacd"
        k = parameters(4);
        theta = psi/gamma(1+1/k);
        z=1-exp(-(xi./theta).^(k));
    elseif model=="bacd"
        kappa=parameters(4);
        sigma2 = parameters(5);
        
        xii = psi.*sigma2^(1+1/kappa).*gamma(1+1/sigma2)./gamma(1+1/kappa)./gamma(1/sigma2-1/kappa);
        z = 1-(1+sigma2.*(xi./xii).^(kappa)).^(-1/sigma2);
    elseif model=="ggacd"
        k = parameters(5);
        delta = parameters(4);
        
        theta = psi.*gamma(k)/gamma(k+1/delta);
        z = gamcdf((xi./theta).^(delta),k,1);
    elseif model=="led1"
        nu0 = 1;
        nu1 = parameters(4);
        [gamma0,gamma1,gamma2,~,~,~,~,~,~] = mygamma(nu0,nu1,0,0,0);
        
        mux = mymu(nu0,nu1,0,0,0);
        z = gamma2/2*myint(xi./psi./mux,2)-(gamma1+2*gamma2)*myint(xi./psi./mux,1)+(gamma0+gamma1+gamma2)*myint(xi./psi./mux,0);
    elseif model =="led2"
        nu0 = 1;
        nu1 = parameters(4);
        nu2 = parameters(5);
        
        mux = mymu(nu0,nu1,nu2,0,0);
        [gamma0,gamma1,gamma2,gamma3,gamma4,~,~,~,~] = mygamma(nu0,nu1,nu2,0,0);
        z = gamma4/24*myint(xi./psi./mux,4)-(2*gamma4/3+gamma3/6)*myint(xi./psi./mux,3)+(gamma2/2+3*gamma3/2+3*gamma4)*...
            myint(xi./psi./mux,2)-(gamma1+2*gamma2+3*gamma3+4*gamma4)*...
            myint(xi./psi./mux,1)+(gamma0+gamma1+gamma2+gamma3+gamma4)*myint(xi./psi./mux,0);
    elseif model == "led3"
        nu0 = 1;
        nu1 = parameters(4);
        nu2 = parameters(5);
        nu3 = parameters(6);
        [gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,~,~] = mygamma(nu0,nu1,nu2,nu3,0);
        
        mux = mymu(nu0,nu1,nu2,nu3,0);
        z = gamma6/720*myint(xi./psi./mux,6)-(gamma6/20+gamma5/120)*myint(xi./psi./mux,5)+(gamma4/24+5*gamma5/24+5*gamma6/8)*myint(xi./psi./mux,4)...
            -(gamma3/6+2*gamma4/3+5*gamma5/3+10*gamma6/3)*myint(xi./psi./mux,3)+(gamma2/2+3*gamma3/2+3*gamma4+5*gamma5+15*gamma6/2)*myint(xi./psi./mux,2)...
            -(gamma1+2*gamma2+3*gamma3+4*gamma4+5*gamma5+6*gamma6)*myint(xi./psi./mux,1)+(gamma0+gamma1+gamma2+gamma3+gamma4+gamma5+gamma6)*myint(xi./psi./mux,0);
    else
        nu0 = 1;
        nu1 = parameters(4);
        nu2 = parameters(5);
        nu3 = parameters(6);
        nu4 = parameters(7);
        [gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8] = mygamma(nu0,nu1,nu2,nu3,nu4);
        
        mux = mymu(nu0,nu1,nu2,nu3,nu4);
        z = gamma8/factorial(8)*myint(xi./psi./mux,8)-(gamma8/630+gamma7/factorial(7))*myint(xi./psi./mux,7)+(gamma6/720+7*gamma7/720+7*gamma8/180)*myint(xi./psi./mux,6)...
            -(7*gamma8/15+7*gamma7/40+gamma6/20+gamma5/120)*myint(xi./psi./mux,5)+(35*gamma8/12+35*gamma7/24+gamma4/24+5*gamma5/24+5*gamma6/8)*myint(xi./psi./mux,4)...
            -(gamma3/6+2*gamma4/3+5*gamma5/3+10*gamma6/3+35*gamma7/6+28*gamma8/3)*myint(xi./psi./mux,3)+(gamma2/2+3*gamma3/2+3*gamma4+5*gamma5+15*gamma6/2+21*gamma7/2+14*gamma8)*myint(xi./psi./mux,2)...
            -(gamma1+2*gamma2+3*gamma3+4*gamma4+5*gamma5+6*gamma6+7*gamma7+8*gamma8)*myint(xi./psi./mux,1)+(gamma0+gamma1+gamma2+gamma3+gamma4+gamma5+gamma6+gamma7+gamma8)*myint(xi./psi./mux,0);
    end
end
function [gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8] = mygamma(nu0,nu1,nu2,nu3,nu4)
    gamma0 = 1;
    gamma1 = (2*nu0*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu4^2)...
        /(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma2 = (2*nu0*nu2-8*nu1*nu2+6*nu1*nu3-24*nu2*nu3+12*nu2*nu4-48*nu3*nu4+2*nu1^2+...
        10*nu2^2+24*nu3^2+44*nu4^2)/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma3 = (2*nu0*nu3+6*nu1*nu2-12*nu1*nu3+8*nu1*nu4+48*nu2*nu3-48*nu2*nu4+156*nu3*nu4-...
        12*nu2^2-56*nu3^2-152*nu4^2)/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma4 = (2*nu0*nu4+8*nu1*nu3-16*nu1*nu4-48*nu2*nu3+88*nu2*nu4-304*nu3*nu4+6*nu2^2+...
        78*nu3^2+346*nu4^2)/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma5 = (10*nu1*nu4+20*nu2*nu3-80*nu2*nu4+360*nu3*nu4-60*nu3^2-520*nu4^2)/...
        (nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma6 = (30*nu2*nu4-240*nu3*nu4+20*nu3^2+500*nu4^2)/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma7 = (70*nu3*nu4-280*nu4^2)/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
    gamma8 = 70*nu4^2/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^4);
end
function y = myint(xi,n)
    res=zeros(length(xi),1);
    for j=0:n
        res = res + factorial(n)/factorial(n-j)*xi.^(n-j).*exp(-xi);
    end
    y = factorial(n)-res;
end
function mux = mymu(nu0,nu1,nu2,nu3,nu4)
    mux = 1/(1-(2*nu0*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu4^2)/(nu0^2+nu1^2+nu2^2+nu3^2+nu4^2));
end

