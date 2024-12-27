function h = weibull_hazard(x,k)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
lambda=1/gamma(1+1/k);
f=k/lambda*(x/lambda).^(k-1).*exp(-(x/lambda).^k);
F=1-exp(-(x/lambda).^k);
h=f./(1-F);
end