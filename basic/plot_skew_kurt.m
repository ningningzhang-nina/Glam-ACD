clear all
load("skew_kurt.mat")
load("new_skew_kurt.mat")
sk1=results1(:,3:4);
[xx,yy1,yy2]=myboundary(sk1((sk1(:,1)<10)&(sk1(:,2)<100),1),sk1((sk1(:,1)<10)&(sk1(:,2)<100),2));
sk11=[[fliplr(xx)';xx'],[fliplr(yy1)';yy2']];
sk2=results2(:,3:4);
[xx,yy1,yy2]=myboundary(sk2((sk2(:,1)<10)&(sk2(:,2)<100),1),sk2((sk2(:,1)<10)&(sk2(:,2)<100),2));
sk22=[[fliplr(xx)';xx'],[fliplr(yy1)';yy2']];
sk3=results3(:,5:6);
sk3=sk3((abs(sk3(:,1))<10)&(sk3(:,2)<100)&(sk3(:,2)>0)&(sk3(:,2)>=sk3(:,1).^2+1),:);
sk31=sk3(sk3(:,1)<0,:);
k = boundary(sk31(:,1),sk31(:,2));
sk331=[sk31(k,1),sk31(k,2)];
sk32=sk3(sk3(:,1)>0,:);
k = boundary(sk32(:,1),sk32(:,2));
sk332=[sk32(k,1),sk32(k,2)];
%%
sk4=results4(:,6:7);
sk4=sk4((abs(sk4(:,1))<10)&(sk4(:,2)<100)&(sk4(:,2)>0)&(sk4(:,2)>=sk4(:,1).^2+1),:);
sk41=sk4(sk4(:,1)<0,:);
k = boundary(sk41(:,1),sk41(:,2));
sk441=[sk41(k,1),sk41(k,2)];
sk42=sk4(sk4(:,1)>0,:);
k = boundary(sk42(:,1),sk42(:,2));
sk442=[sk42(k,1),sk42(k,2)];
%%
new441=[];
for i=1:length(sk441)
    if sk441(i,2)<sk441(i,1)*(-100/9)
        new441=[new441;sk441(i,:)];
    end
end
new442=[];
for i=1:length(sk442)
    if sk442(i,2)<sk442(i,1)*(100/9)
        new442=[new442;sk442(i,:)];
    end
end
%%
figure('Renderer', 'painters', 'Position', [10 10 900 600],'DefaultAxesFontSize',18)
plot(sk11(:,end-1),sk11(:,end),LineWidth=1.5,color=[0.35 0.35 0.35])
hold on
% fill(sk11(:,end-1),sk11(:,end),'r')
% hold on
plot(sk22(:,end-1),sk22(:,end),'-.',LineWidth=1.5,color=[0.35 0.35 0.35])
% hold on
% fill(sk22(:,end-1),sk22(:,end),'g')
hold on
plot(sk331(:,end-1),sk331(:,end),':',LineWidth=1.5,color=[0.35 0.35 0.35])
hold on

plot(new441(:,end-1),new441(:,end),'--',LineWidth=1.5,color=[0.35 0.35 0.35])

xxx=-10:0.1:10;
yyy=xxx.^2+1;
hold on
plot(xxx,yyy,'-o',color=[0.35 0.35 0.35])
hold on
plot(sk332(:,end-1),sk332(:,end),':',LineWidth=1.5,color=[0.35 0.35 0.35])
hold on
plot(new442(:,end-1),new442(:,end),'--',LineWidth=1.5,color=[0.35 0.35 0.35])
hold on
legend('Generalized gamma','Burr','Glam_1','Glam_2','Skewness-Kurtosis frontier','Location','southoutside','Orientation','horizontal')
ylim([0,55])
xlim([-6,6])
xlabel("Skewness")
ylabel("Kurtosis")
saveas(gcf,'skew_kurt_compare','epsc')
