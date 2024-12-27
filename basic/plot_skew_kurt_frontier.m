clear all
x=1;
%%
results1=[];
a=1;
for delta=[0.1:0.01:3,3:1:100]
    for kappa=0.1:0.5:100
        [~,~,~,s,k]=compute_moments(x,"ggacd",[a,delta,kappa]);
        results1=[results1;[delta,kappa,s,k]];
    end
end
sk1=results1(:,3:4);
[xx,yy1,yy2]=myboundary(sk1((sk1(:,1)<15)&(sk1(:,2)<300),1),sk1((sk1(:,1)<15)&(sk1(:,2)<300),2));
sk11=[[fliplr(xx)';xx'],[fliplr(yy1)';yy2']];
%%
results2=[];
lambda=1;

for sigma2=0.01:0.01:1
    for kappa=sigma2*4+0.01:0.1:100
        [~,~,~,s,k]=compute_moments(x,"bacd",[lambda,kappa,sigma2]);
        results2=[results2;[sigma2,kappa,s,k]];
    end
end
sk2=results2(:,3:4);
[xx,yy1,yy2]=myboundary(sk2((sk2(:,1)<4)&(sk2(:,2)<160),1),sk2((sk2(:,1)<4)&(sk2(:,2)<160),2));
sk22=[[fliplr(xx)';xx'],[fliplr(yy1)';yy2']];
%%
x=1;
k=1;
new_results3=[];
for nu1=[-20:1:0,0:0.1:5,5:1:20]
    tic;
    for s=[0.1:0.1:10,10:10:50]
        for alpha=[-0.9:0.01:0.1,0.1:10:50]
            [f,F,h,moments]=new_gled_prop([alpha,s,k,nu1],x);
            new_results3=[new_results3;alpha,s,k,nu1,moments(3),moments(4)];
        end
    end
    toc;
end
%%
x=1;
k=1;
%results4=[];
for nu2=[0,5,10]
    tic;
    for nu1=[-20:1:0,0:0.1:5,5:1:20]
        for s=[0.1:0.1:10,10:5:100]
            for alpha=[-0.9:0.01:0.1,0.1:5:100]
                [f,F,h,moments]=new_gled_prop([alpha,s,k,nu1,nu2],x);
                results4=[results4;alpha,s,k,nu1,nu2,moments(3),moments(4)];
            end
        end
    end
    toc;
end
%%
figure
plot(results4(:,end-1),results4(:,end),'o')
hold on
plot(results3(:,end-1),results3(:,end),'o')
hold on
plot(results2(:,end-1),results2(:,end),'o')
hold on
plot(results1(:,end-1),results1(:,end),'o')
xlim([-12,12])
ylim([0,100])
%%
figure
plot(results0(:,end-1),results0(:,end),'o')
%%
x=1;
k=1;
results1=[];
tic;
for alpha=[-0.9,0,1,2,20]
    for s=[0.5,10]
        for nu1=-10:0.01:10
            [f,F,h,moments]=new_gled_prop([alpha,s,k,nu1],x);
            results1=[results1;alpha,s,k,nu1,moments(3),moments(4)];
        end
    end
end
toc;


%%
figure
plot(results1(:,end-1),results1(:,end),'o')
hold on
plot(results2(:,end-1),results2(:,end),'o')
hold on
plot(results3(:,end-1),results3(:,end),'o')
hold on
plot(results4(:,end-1),results4(:,end),'o')
xlim([-10,10])
ylim([0,100])