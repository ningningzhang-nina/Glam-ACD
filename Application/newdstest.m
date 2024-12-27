function [pt,f,fpar,pvalues]=newdstest(xi,xi1,psi1,parameters,model,type)
    %[h,f,~,~]=kde(xi,length(xi),[],[]);
    pt=log(xi);
    [f,~,h]=ksdensity(pt,pt);
    fpar = newcompute_density(xi,xi1,psi1,parameters,model,type);
    %fpar=xi.*fpar;

    expDs = mean(f)/2/sqrt(pi*h);
    varDs = mean(f.^3)/2/sqrt(2*pi);
    DS = length(xi)*sqrt(h)*mean((fpar-f).^2);

    pvalues = 2*(1-cdf('Normal',abs(DS-expDs)/sqrt(varDs),0,1));
end
function fpar = newcompute_density(xi,xi1,psi1,parameters,model,type)
    if type==1
        if (model=="eacd")||(model=="logeacd1")||(model=="logeacd2")
            fpar = xi./psi1.*exp(-xi./psi1);
        elseif (model=="wacd")||(model=="logwacd1")||(model=="logwacd2")
            kappa = parameters(4);

            fpar = xi./psi1.*kappa.*((xi./psi1.*gamma(1+1/kappa)).^(kappa-1)).*exp(-(xi./psi1.*gamma(1+1/kappa)).^(kappa));
        elseif (model=="ggacd")||(model=="logggacd1")||(model=="logggacd2")
        
            k = parameters(5);
            delta = parameters(4);
            thetai = psi1.*gamma(k)/gamma(k+1/delta);
            fpar = xi./thetai.*delta.*(xi.^(delta*k-1))./(thetai.^(delta*k-1))./gamma(k).*exp(-(xi./thetai).^delta);
        elseif (model=="bacd")||(model=="logbacd1")||(model=="logbacd2")
       
            kappa=parameters(4);
            sigma2 = parameters(5);
            xii = psi1.*sigma2.^(1+1/kappa)*gamma(1/sigma2+1)/gamma(1+1/kappa)/gamma(1/sigma2-1/kappa);
            fpar = xi./psi1.*kappa.*(xii.^(1-kappa)).*(xi.^(kappa-1))./((1+sigma2.*(xi./xii).^(kappa)).^(1/sigma2+1));
        elseif (model=="led1")||(model=="logled11")||(model=="logled12")
            nu0=1;
            nu1=parameters(4);
            fpar = led(xi,psi1,nu0,nu1,0,0,0);
        elseif (model=="led2")||(model=="logled21")||(model=="logled22")
            nu0=1;
            nu1=parameters(4);
            nu2=parameters(5);
            fpar = led(xi,psi1,nu0,nu1,nu2,0,0);
        elseif (model=="led3")||(model=="logled31")||(model=="logled32")
            nu0=1;
            nu1=parameters(4);
            nu2=parameters(5);
            nu3=parameters(6);
            fpar = led(xi,psi1,nu0,nu1,nu2,nu3,0);
        elseif (model=="led4")||(model=="logled41")||(model=="logled42")
            nu0=1;
            nu1=parameters(4);
            nu2=parameters(5);
            nu3=parameters(6);
            nu4=parameters(7);
            fpar = led(xi,psi1,nu0,nu1,nu2,nu3,nu4);
        end
    else
        alpha = parameters(1);
        beta = parameters(2);
        Gamma = parameters(3);
        if model=="eacd"
            psi=zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1=psi(i);
            end
        fpar = xi./psi.*exp(-xi./psi);
        elseif model=="logeacd1"
            psi=zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*log(xi1(i))+Gamma*log(psi1);
                psi1=exp(psi(i));
            end
        fpar = xi./exp(psi).*exp(-xi./exp(psi));
        elseif model=="logeacd2"
            psi=zeros(length(xi),1);
            epsilon1=xi1(1)./psi1;
            for i=1:length(xi)
                psi(i) = alpha + beta*epsilon1+Gamma*log(psi1);
                epsilon1=xi(i)/exp(psi(i));
                psi1 = exp(psi(i));
            end
        fpar = xi./exp(psi).*exp(-xi./exp(psi));
        elseif model=="wacd"
            kappa = parameters(4);
            theta = zeros(length(xi),1);
            for i=1:length(xi)
                theta(i) = alpha/gamma(1+1/kappa)+beta/gamma(1+1/kappa)*xi1(i)+Gamma*psi1/gamma(1+1/kappa);
                psi1 = theta(i)*gamma(1+1/kappa);
            end
            fpar = xi./theta.*kappa.*((xi./theta).^(kappa-1)).*exp(-(xi./theta).^(kappa));
        elseif model=="logwacd1"
            kappa = parameters(4);
            theta = zeros(length(xi),1);
            for i=1:length(xi)
                theta(i) = alpha/gamma(1+1/kappa)+beta/gamma(1+1/kappa)*log(xi1(i))+Gamma*log(psi1)/gamma(1+1/kappa);
                psi1 = exp(theta(i))*gamma(1+1/kappa);
            end
            fpar = xi./exp(theta)./gamma(1+1/kappa).*kappa.*((xi./exp(theta)).^(kappa-1)).*exp(-(xi./exp(theta)).^(kappa));
        elseif model=="logwacd2"
            kappa = parameters(4);
            theta = zeros(length(xi),1);
            epsilon1 = xi1(1)./psi1;
            for i=1:length(xi)
                theta(i) = alpha/gamma(1+1/kappa)+beta/gamma(1+1/kappa)*epsilon1+Gamma*log(psi1)/gamma(1+1/kappa);
                epsilon1 = xi(i)/(exp(theta(i))*gamma(1+1/kappa));
                psi1 = exp(theta(i))*gamma(1+1/kappa);
            end
            fpar = xi./exp(theta)./gamma(1+1/kappa).*kappa.*((xi./exp(theta)).^(kappa-1)).*exp(-(xi./exp(theta)).^(kappa));
        elseif model=="bacd"
            kappa=parameters(4);
            sigma2 = parameters(5);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1=psi(i);
            end
            xii = psi.*sigma2^(1+1/kappa).*gamma(1+1/sigma2)./gamma(1+1/kappa)./gamma(1/sigma2-1/kappa);
            fpar = xi./xii.*kappa.*(xii.^(1-kappa)).*(xi.^(kappa-1))./((1+sigma2.*(xi./xii).^(kappa)).^(1/sigma2+1));
        elseif model=="logbacd1"
            kappa=parameters(4);
            sigma2 = parameters(5);
            xii = zeros(length(xi),1);
            korr = sigma2^(1+1/kappa)*gamma(1+1/sigma2)/gamma(1+1/kappa)/gamma(1/sigma2-1/kappa);
            xii1 = log(psi1*korr);
            for i=1:length(xi)
                xii(i) = alpha*korr+beta*korr*log(xi1(i))+Gamma*xii1;
                xii1 = xii(i);
            end
            %xii = exp(psi).*sigma2^(1+1/kappa).*gamma(1+1/sigma2)./gamma(1+1/kappa)./gamma(1/sigma2-1/kappa);
            fpar = xi./exp(xii).*korr.*kappa.*(exp(xii).^(1-kappa)).*(xi.^(kappa-1))./((1+sigma2.*(xi./exp(xii)).^(kappa)).^(1/sigma2+1));
        elseif model=="logbacd2"
            kappa=parameters(4);
            sigma2 = parameters(5);
            xii = zeros(length(xi),1);
            korr = sigma2^(1+1/kappa)*gamma(1+1/sigma2)/gamma(1+1/kappa)/gamma(1/sigma2-1/kappa);
            xii1 = log(psi1*korr);
            
            epsilon1 = xi1(1)./psi1;
            for i=1:length(xi)
                xii(i) = alpha*korr+beta*korr*epsilon1+Gamma*xii1;
                xii1 = xii(i);
                epsilon1 = xi(i)/exp(xii(i))*korr;
            end
            %xii = exp(psi).*sigma2^(1+1/kappa).*gamma(1+1/sigma2)./gamma(1+1/kappa)./gamma(1/sigma2-1/kappa);
            fpar = xi./exp(xii).*korr.*(kappa.*exp(xii).^(1-kappa)).*(xi.^(kappa-1))./((1+sigma2.*(xi./exp(xii)).^(kappa)).^(1/sigma2+1));
        elseif model=="ggacd"
            k = parameters(5);
            delta = parameters(4);
            psi=zeros(length(xi),1);
            theta=zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1=psi(i);
                theta(i) = psi(i).*gamma(k)/gamma(k+1/delta);
            end
            fpar = xi./theta.*delta.*(xi.^(delta*k-1))./(theta.^(delta*k-1))/gamma(k).*exp(-(xi./theta).^delta);
        elseif model=="logggacd1"
            k = parameters(5);
            delta = parameters(4);
            theta=zeros(length(xi),1);
            korr = gamma(k)/gamma(k+1/delta);
            theta1=log(psi1*korr);
            for i=1:length(xi)
                theta(i) = alpha*korr+beta*korr*log(xi1(i))+Gamma*theta1;
                theta1=theta(i);
            end
            fpar = xi./exp(theta).*korr.*delta.*(xi.^(delta*k-1))/(exp(theta).^(delta*k-1))/gamma(k)*exp(-(xi./exp(theta)).^delta);
        elseif model=="logggacd2"
            k = parameters(5);
            delta = parameters(4);
            theta=zeros(length(xi),1);
            korr = gamma(k)/gamma(k+1/delta);
            epsilon1 = xi1(1)./psi1;
            theta1=log(psi1*korr);
            z=zeros(length(xi),1);
            for i=1:length(xi)
                theta(i) = alpha*korr+beta*korr*epsilon1+Gamma*theta1;
                theta1=theta(i);
                epsilon1 = xi(i)/exp(theta(i))*korr;
                %theta = exp(psi(i)).*gamma(k)/gamma(k+1/delta);
                z(i) = gamcdf((xi(i)./exp(theta(i))).^(delta),k,1);
            end
            fpar = xi./exp(theta).*korr.*delta.*(xi.^(delta*k-1))/(exp(theta).^(delta*k-1))/gamma(k)*exp(-(xi./exp(theta)).^delta);
        elseif model=="led1"
            nu0 = 1;
            nu1 = parameters(4);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1 = psi(i);
            end
            fpar = led(xi,psi,nu0,nu1,0,0,0);
        elseif model=="logled11"
            nu0 = 1;
            nu1 = parameters(4);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*log(xi1(i))+Gamma*log(psi1);
                psi1 = exp(psi(i));
            end
            fpar = led(xi,exp(psi),nu0,nu1,0,0,0);
        elseif model=="logled12"
            nu0 = 1;
            nu1 = parameters(4);
            psi = zeros(length(xi),1);
            epsilon1 = xi1(1)./psi1;
            for i=1:length(xi)
                psi(i) = alpha+beta*epsilon1+Gamma*log(psi1);
                psi1 = exp(psi(i));
                epsilon1 = xi1(i)/exp(psi(i));
            end
            fpar = led(xi,exp(psi),nu0,nu1,0,0,0);
        elseif model =="led2"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1 = psi(i);
            end
            fpar = led(xi,psi1,nu0,nu1,nu2,0,0);
        elseif model =="logled21"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*log(xi1(i))+Gamma*log(psi1);
                psi1 = exp(psi(i));
            end
            fpar = led(xi,exp(psi1),nu0,nu1,nu2,0,0);
        elseif model =="logled22"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            psi = zeros(length(xi),1);
            epsilon1 = xi1(1)./psi1;
            for i=1:length(xi)
                psi(i) = alpha+beta*epsilon1+Gamma*log(psi1);
                psi1 = exp(psi(i));
                epsilon1 = xi1(i)/exp(psi(i));
            end
            fpar = led(xi,exp(psi1),nu0,nu1,nu2,0,0);
        elseif model == "led3"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            nu3 = parameters(6);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1 = psi(i);
            end
            fpar = led(xi,psi1,nu0,nu1,nu2,nu3,0);
        elseif model == "logled31"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            nu3 = parameters(6);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*log(xi1(i))+Gamma*log(psi1);
                psi1 = exp(psi(i));
            end
            fpar = led(xi,exp(psi1),nu0,nu1,nu2,nu3,0);
        elseif model == "logled32"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            nu3 = parameters(6);
            psi = zeros(length(xi),1);
            epsilon1 = xi1(1)./psi1;
            for i=1:length(xi)
                psi(i) = alpha+beta*epsilon1+Gamma*log(psi1);
                psi1 = exp(psi(i));
                epsilon1 = xi1(i)/exp(psi(i));
            end
            fpar = led(xi,exp(psi1),nu0,nu1,nu2,nu3,0);
        elseif model=="led4"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            nu3 = parameters(6);
            nu4 = parameters(7);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*xi1(i)+Gamma*psi1;
                psi1 = psi(i);
            end
            fpar = led(xi,psi1,nu0,nu1,nu2,nu3,nu4);
        elseif model=="logled41"
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            nu3 = parameters(6);
            nu4 = parameters(7);
            psi = zeros(length(xi),1);
            for i=1:length(xi)
                psi(i) = alpha+beta*log(xi1(i))+Gamma*log(psi1);
                psi1 = exp(psi(i));
            end
            fpar = led(xi,exp(psi1),nu0,nu1,nu2,nu3,nu4);
        else
            
            nu0 = 1;
            nu1 = parameters(4);
            nu2 = parameters(5);
            nu3 = parameters(6);
            nu4 = parameters(7);
            psi = zeros(length(xi),1);
            epsilon1 = xi1(1)./psi1;
            for i=1:length(xi)
                psi(i) = alpha+beta*epsilon1+Gamma*log(psi1);
                psi1 = exp(psi(i));
                epsilon1 = xi1(i)/exp(psi(i));
            end
            fpar = led(xi,exp(psi1),nu0,nu1,nu2,nu3,nu4);
        end
    end
end
function fpar = led(xi,mux,nu0,nu1,nu2,nu3,nu4)
    
    nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2+nu4^2;
    a_nu_i = 1/(1-(2*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu2^2)/nu_i_square);
    fpar = xi.*exp(-xi./mux./a_nu_i).*(nu0*lagu(xi./mux./a_nu_i,0)+nu1*lagu(xi./mux./a_nu_i,0)+nu2*lagu(xi./mux./a_nu_i,0)+...
        nu3*lagu(xi./mux./a_nu_i,0)+nu4*lagu(xi./mux./a_nu_i,0))./mux./a_nu_i;

end
function l= lagu(x,r)
    if r==0
        l=1;
    elseif r==1
        l=-x+1;
    elseif r==2
        l=(x.^2-4*x+1)./2;
    elseif r==3
        l=(-x.^3+9.*x.^2-18.*x+6)./6;
    elseif r==4
        l=(x.^4-16.*x.^3+72.*x^2-96.*x+24)./24;
    end
end