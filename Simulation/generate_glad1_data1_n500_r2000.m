clear all
syms x f(x) z(x)
% initiate parameters
replication_time=2000;
num=500;
alpha=0;
s=1;
nu=-0.8240;
k=(1+nu^2)/(1-2*nu+2*nu^2);
% first parameter setting: mu=1,sigma=2,nu1=-2,nu2=nu3=nu4=0
m=1;
% density and cdf
f(x)=0;
for kk=0:2*m
    f(x)=f(x)+(s/k)*(x./k).^(s-1+s*alpha).*exp(-(x./k).^s)./...
        mystandardlized(nu,alpha).*gammak(nu,alpha,kk).*GL(alpha,kk,(x./k).^s);
end

randx = zeros(num,replication_time);
% generate random number
%%
count=0;
for repeat=1:replication_time
    tic;
    ite=0;
    while ite<=num
        count=count+1;
        try
            rng(count)
            r = random('Uniform',0,1);
            z(x)=0;
            for kk=0:2*m
                for mm=0:kk
                    z(x)=z(x)+1/mystandardlized(nu,alpha)*gammak(nu,alpha,kk)*(-1)^(mm)*...
                    mynchoosek(kk+alpha,kk-mm)/factorial(mm)*(gamma(alpha+mm+1)-igamma(alpha+mm+1,(x.^2./k).^s));
                end
            end
            z(x)=z(x)-r;
            a=(vpasolve(z==0,x)).^2;
            ite=ite+1;
            randx(ite,repeat) =a;
            
        catch
            warning("not solved")
        end
    end

    % [parameters] = my_glgd(randx(1:end,repeat),[alpha,s,k,nu]);
    % final_results{}(repeat,:)=parameters;
    toc;
    % repeat 
    % final_results(repeat,:)-[alpha,s,k,nu]
end

%mean(final_results(:,1:4)-[alpha,s,k,nu])
save("simulation_glgd1_data1_n500_r2000.mat")
function a=mynchoosek(n,r)
    a=1;
    for i=1:r
        a=a*(n+1-i)/i;
    end
    a=a;
end

function L=GL(alpha,n,x)
    L=0;
    for m=0:n
        L=L+(-1)^m*mynchoosek(n+alpha,n-m).*x.^m./factorial(m);
    end
end
% function I=myintegral(m,alpha,x)
% I=gamma(m+alpha+1);
% for j=0:(m+alpha)
%     I=I-gamma(m+alpha+1)/gamma(m+alpha-j+1)*x.^(m+alpha-j).*exp(-x);
% end
% end
function g=gammak(nu,alpha,l)
m=length(nu);
aijl=zeros(m+1,m+1);

for i=0:m
    for j=0:m
        if (l>=abs(i-j))&&(l<=i+j)
            c=0;
            for r=max(max(i,j),l):floor((i+j+l)/2)
                c=c+gamma(alpha+r+1)*(-2)^(i+j+l-2*r)/factorial(r-i)/factorial(r-j)/...
                    factorial(r-l)/factorial(i+j+l-2*r);
            end
            aijl(i+1,j+1)=factorial(l)/gamma(alpha+l+1)*c;
        else
            aijl(i+1,j+1)=0;
        end
    end
end
g=[1,nu]*aijl*[1,nu]';
end
function s=mystandardlized(nu,alpha)
m=length(nu);
s=gamma(1+alpha);
for i=1:m
s=s+nu(i)^2*gamma(alpha+1+i)/factorial(i);
end
end
function mu=compute_noncentral_moments(nu,alpha,n,s,k)
m=length(nu);
mu=0;
for kk=0:2*m
    for mm=0:kk
        mu=mu+k^n/mystandardlized(nu,alpha)*gammak(nu,alpha,kk)*(-1)^(mm)*...
                mynchoosek(kk+alpha,kk-mm)/factorial(mm)*gamma(alpha+mm+1+n/s);
    end
end
end
function [hx,sx,kx]=compute_high_moments(mu1,mu2,mu3,mu4)
hx = -mu1.^2+mu2;
sx=(2*mu1.^3-3*mu1*mu2+mu3)/(hx^(3/2));
kx=(-3*mu1.^4+6*mu1.^2.*mu2-4*mu1.*mu3+mu4)/(hx^(2));
end