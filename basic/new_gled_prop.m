function [f,F,h,moments]=new_gled_prop(params,x)
    alpha=params(1);
    k=params(3);
    s=params(2);
    nu=params(4:end);
    % s 控制x的范围，s越小，x的范围越大
    % k也是控制x的范围，k越大，x的范围越大
    m=length(nu);
    f=0;
    for kk=0:2*m
        f=f+(s/k)*(x./k).^(s-1+s*alpha).*exp(-(x./k).^s)./...
            mystandardlized(nu,alpha).*gammak(nu,alpha,kk).*GL(alpha,kk,(x./k).^s);
    end
    
    %f=fun(alpha,s,k,x);
    F=0;
    for kk=0:2*m
        for mm=0:kk
            %fun=@(xx)(s/k)*(xx./k).^(s-1+s*alpha+s*m).*exp(-(xx./k).^s);
            F=F+1/mystandardlized(nu,alpha)*gammak(nu,alpha,kk)*(-1)^(mm)*...
                mynchoosek(kk+alpha,kk-mm)/factorial(mm)*gammainc((x./k).^s,alpha+mm+1)*gamma(alpha+mm+1);
            %myintegral(mm,alpha,(x./k).^s);
        end
    end
    h=f./(1-F);
    mu=compute_noncentral_moments(nu,alpha,1,s,k);
    m2=compute_noncentral_moments(nu,alpha,2,s,k);
    m3=compute_noncentral_moments(nu,alpha,3,s,k);
    m4=compute_noncentral_moments(nu,alpha,4,s,k);
    [hx,sx,kx]=compute_high_moments(mu,m2,m3,m4);
    if (mu>0)&&(hx>0)&&(kx>0)&&(kx-sx^2-1>0)
        moments=[mu;hx;sx;kx];
    else
        moments=[Inf;Inf;Inf;Inf];
    end
end
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
muu=0;
for kk=0:2*m
    for mm=0:kk
        muu=muu+gammak(nu,alpha,kk)*(-1)^(mm)*...
                mynchoosek(kk+alpha,kk-mm)/factorial(mm)*gamma(alpha+mm+1+n/s);
    end
end
mu=k^n/mystandardlized(nu,alpha)*muu;
end
function [hx,sx,kx]=compute_high_moments(mu1,mu2,mu3,mu4)
hx = -mu1.^2+mu2;
sx=(2*mu1.^3-3*mu1*mu2+mu3)/(hx^(3/2));
kx=(-3*mu1.^4+6*mu1.^2.*mu2-4*mu1.*mu3+mu4)/(hx^(2));
end