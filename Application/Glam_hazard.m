function h = Glam_hazard(x,nu)
    alpha=0;
    s=1;
    
    m=length(nu);
    if m==1
        nu_i_square = 1+nu(1)^2;
        k = 1/(1-(2*nu(1)-2*nu(1)^2)/nu_i_square);
    elseif m==2
        nu_i_square = 1+nu(1)^2+nu(2)^2;
        k = 1/(1-(2*nu(1)+4*nu(1)*nu(2)-2*nu(1)^2-4*nu(2)^2)/nu_i_square);
    elseif m==3
        nu_i_square = 1+nu(1)^2+nu(2)^2+nu(3)^2;
        k = 1/(1-(2*nu(1)+4*nu(1)*nu(2)+6*nu(2)*nu(3)-2*nu(1)^2-4*nu(2)^2-6*nu(3)^2)/nu_i_square);
    else
        nu_i_square = 1+nu(1)^2+nu(2)^2+nu(3)^2+nu(4)^2;
        k=1/(1-(2*nu(1)+4*nu(1)*nu(2)+6*nu(2)*nu(3)+8*nu(3)*nu(4)-2*nu(1)^2-4*nu(2)^2-6*nu(3)^2-8*nu(2)^2)/nu_i_square);

    % else
    %     k=(1+sum(nu.^2))/(1-gamma1(1,nu(1),nu(2),nu(3),nu(4)));
    end
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



% function gamma_1 = gamma1(nu0,nu1,nu2,nu3,nu4)
% 
% nu_i_square = nu0^2+nu1^2+nu2^2+nu3^2+nu4^2;
% gamma_1 = (2*nu0*nu1+4*nu1*nu2+6*nu2*nu3+8*nu3*nu4-2*nu1^2-4*nu2^2-6*nu3^2-8*nu4^2)/nu_i_square;
% end