function [Lower_bound,Upper_bound]=generate_confidence_interval(z)
    f=length(z);
    conf_interval=0.95;
    bindwith=1/20;
    m=1:1:f;
    p=zeros(f+1,1);
    q=1;
    p(1)=(1-bindwith)^f;
    p(f+10)=(bindwith)^f;
    while q <= f
      p(q+1,1)=exp(sum(log(m(1:length(m))))-sum(log(m(1:length(m)-q)))-sum(log(m(1:q)))+ q*log(bindwith)+(f-q)*log(1-bindwith));
      q=q+1;
    end
    
    s=1;
    kum=p;
    while s <=f
    kum(s+1,1)=kum(s,1)+p(s+1,1);
    s=s+1;
    end
    t=1;
    g=0;
    while t <= f+1
    if kum(t,1) <= (1-conf_interval)/2
        g=g+1;
    else
        g=g+0;
    end
    t=t+1;
    end
    Lower_bound=(g-1)/f;
    r=1;
    v=0;
    while r <= f+1
    if kum(r,1) <= (1-(1-conf_interval)/2)
      v=v+1;
    else
      v=v+0;
    end
      r=r+1;
    end
    Upper_bound=v/f;
end
