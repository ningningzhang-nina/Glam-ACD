function ft = gamma_ksdensity(xi,pts,b)
%x: sample data
%pts: points at which to evaluate f
n = length(xi);
ft = zeros(length(pts),1);

for i = 1:length(pts)
    ft(i) = sum(gamma_kernel2(pts(i),b,xi))/n;
end
%     function k = gamma_kernel1(x,b,t)
%         k = t.^(x/b).*exp(-t./b)./b.^(x./b+1)./gamma(x./b+1);
%     end
    function k = gamma_kernel2(x,b,t)
        if x>=2*b
            ro = x./b;
        else
            ro = (x./b).^2./4+1;
        end
        k = t.^(ro-1).*exp(-t./b)./b.^(ro)./gamma(ro);
    end
end