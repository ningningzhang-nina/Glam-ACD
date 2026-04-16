function [c,p]=goodnesstest(z,counts)
    n = length(z);
%     pi = zeros(length(counts),1);
%     for i=1:length(counts)
% 
%         pi(i) = nchoosek(n,counts(i))*(1/20)^(counts(i))*(1-1/20)^(n-counts(i));
%     end
    c=0;
    for j=1:20
        c=c+(counts(j)-n/20)^2/(n/20);
    end
    p = 1-chi2cdf(c,19);
end