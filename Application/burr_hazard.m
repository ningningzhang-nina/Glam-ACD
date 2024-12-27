function h = burr_hazard(x,parameters)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
kappa=parameters(1);
sigma2=parameters(2);
h=kappa*x.^(kappa-1)./(1+sigma2*x.^(kappa));
end