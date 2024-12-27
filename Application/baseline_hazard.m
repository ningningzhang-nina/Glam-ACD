function [t,h] = baseline_hazard(x,k)
% input duration series: x1,x2...xn
t=sort(x);
%h=zeros(1,length(t));
i=k+1;
while i<=length(t)-k
    h(i)=2*k/((length(x)-i)*(t(i+k)-t(i-k)));
    i=i+1;
end
t=t(k+1:i-1);
h=h(k+1:end);
end