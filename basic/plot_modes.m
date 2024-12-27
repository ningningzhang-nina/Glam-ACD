clear all
figure('Renderer','painters','Position',[10,10,1500,900],'DefaultAxesFontSize',22)
%% no mode
subplot(5,2,1)
x=0:0.1:10;
alpha=0;
s=1;
k=1;
nu1=0;
params=[alpha,s,k,nu1];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f,'DisplayName',"\boldmath$\theta$=(0,1,1,0,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("No mode")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% one mode
subplot(5,2,2)
x=0:0.1:10;
alpha=0;
s=2;
k=2;
nu1=-1;
params=[alpha,s,k,nu1];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,2,2,$-$1,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("One mode")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% one mode
subplot(5,2,3)
x=0:0.001:35;
alpha=0;
s=0.7966;
k=2;
nu1=-0.6467;
params=[alpha,s,k,nu1];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,2,0.7966,$-$0.6467,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("One mode")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% two modes
subplot(5,2,4)
x=0:0.1:10;
alpha=0;
s=2;
k=2;
nu1=1;
params=[alpha,s,k,nu1];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,2,2,1,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Two modes")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% two modes
subplot(5,2,5)
x=0:0.1:50;
alpha=0;
s=1.8027;
k=5;
nu1=-0.6917;
nu2=0.4801;
nu3=-0.4734;
nu4=0.0675;
params=[alpha,s,k,nu1,nu2,nu3,nu4];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,5,1.8027,$-$0.6917,0.4801,$-$0.4734,0.0675)",LineWidth=2,color=[0.35 0.35 0.35])
ylim([0,0.10])
title("Two modes")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% three modes
subplot(5,2,6)
x=0:0.1:10;
alpha=0;
s=2;
k=2;
nu1=1;
nu2=10;
params=[alpha,s,k,nu1,nu2];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,2,2,1,10,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Three modes")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% three modes
subplot(5,2,7)
x=0:0.1:30;
alpha=0;
s=2.5937;
k=5;
nu1=-0.4075;
nu2=0.6181;
nu3=-0.2995;
nu4=0.5788;
params=[alpha,s,k,nu1,nu2,nu3,nu4];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,5,2.5937,$-$0.4075,0.6181,$-$0.2995,0.5788)",LineWidth=2,color=[0.35 0.35 0.35])
ylim([0,0.18])
title("Three modes")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% four modes
subplot(5,2,8)
x=0:0.1:20;
alpha=0;
s=1;
k=1;
nu1=2;
nu2=1;
nu3=2;
nu4=-5;
params=[alpha,s,k,nu1,nu2,nu3,nu4];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,1,1,2,1,2,$-$5)",LineWidth=2,color=[0.35 0.35 0.35])
title("Four modes")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% four modes
subplot(5,2,9)
x=0.01:0.01:80;
alpha=0;
s=1.1606;
k=5;
nu1=0.2243;
nu2=-0.6619;
nu3=-0.1279;
nu4=0.5385;
params=[alpha,s,k,nu1,nu2,nu3,nu4];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(0,5,1.1606,0.2243,$-$0.6619,$-$0.1279,0.5385)",LineWidth=2,color=[0.35 0.35 0.35])
ylim([0,0.08])
title("Four modes")
hl = legend('show');
set(hl, 'Interpreter','latex')
hold on
%% five modes
subplot(5,2,10)
x=0:0.1:15;
alpha=1;
s=2;
k=2;
nu1=2;
nu2=-1;
nu3=-2;
nu4=-2;
params=[alpha,s,k,nu1,nu2,nu3,nu4];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,f, 'DisplayName',"\boldmath$\theta$=(1,2,2,2,$-$1,$-$2,$-$2)",LineWidth=2,color=[0.35 0.35 0.35])
title("Five modes")%,'FontWeight', 'bold'
hl = legend('show');
set(hl, 'Interpreter','latex')
saveas(gcf,"modes.eps")
