function [xx,low,upper]=myboundary(x,y)
    xx=min(x):0.2:max(x);
    for i=1:length(xx)
        try
        k=find((x<=xx(i)+0.3)&(x>=xx(i)-0.3));
        catch
            k=find((x<=xx(i)+2)&(x>=xx(i)-2));
        end
        low(i)=min(y(k));
       
        upper(i)=max(y(k));
    end

end
