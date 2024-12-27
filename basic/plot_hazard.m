clear all
figure('Renderer','painters','Position',[10,10,1500,900],'DefaultAxesFontSize',22)
%% 
subplot(5,2,1)
x=0:0.01:10;
alpha=0;
s=1;
k=1;
nu1=0;
params=[alpha,s,k,nu1];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h,'DisplayName',"\boldmath$\theta$=(0,1,1,0,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
ylim([0.9,1.1])
title("Constant")
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')
hold on
%% increasing
subplot(5,2,2)
x=0:0.01:3;
% alpha=0;
% s=2;
% k=2;
% nu1=-1;
params=[0.3446,0.9338,1,0,0,0,0];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=(0.3446,0.9338,1,0,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Increasing")
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')
hold on
%% decreasing
subplot(5,2,3)
x=0:0.01:3;
% alpha=1;
% s=1/2;
% k=1;
% nu1=0;
params=[0.0868,0.3449,1,0,0,0,0];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=(0.0868,0.3449,1,0,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Decreasing")
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')
hold on

%% U-shaped
subplot(5,2,4)
x=0:0.01:3;
% alpha=0;
% s=1;
% k=1;
% nu1=2;
params=[-0.5948,1.4934,1,0,0,0,0];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=($-$0.5948,1.4934,1,0,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("U-shaped")
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')
hold on
%% Inverted U-shaped
subplot(5,2,5)
x=0:0.01:20;
% alpha=1;
% s=0.8;
% k=1;
% nu1=0;
params=[0.5426,0.7885,1,0,0,0,0];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=(0.5426,0.7885,1,0,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Inverted U-shaped")
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')
hold on
% One mode
subplot(5,2,6)
x=0:0.01:7;
% alpha=0;
% s=2;
% k=2;
% nu1=1;
params=[-0.6952,1.6516,1,0.0767];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=($-$0.6952,1.6516,1,0.0767,0,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
ylim([0,15])
title("Three extrema")
hl = legend('show','Location','northwest');
set(hl, 'Interpreter','latex')
hold on
%  %
% subplot(5,2,6)
% x=0:0.01:20;
% alpha=0;
% s=1;
% k=1;
% nu1=2;
% nu2=-10;
% params=[alpha,s,k,nu1,nu2];
% [f,F,h,moments]=new_gled_prop(params,x);
% 
% plot(x,h, 'DisplayName',"\boldmath$\theta$=(0,1,1,2,$-$10,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
% title("Three extrema")
% hl = legend('show','Location','southeast');
% set(hl, 'Interpreter','latex')
% hold on
%%
subplot(5,2,7)
x=0:0.01:7;
% alpha=0;
% s=2;
% k=2;
% nu1=1;
% nu2=10;
params=[-0.9245,1.7703,1,0.8266,0.5924,0,0];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=($-$0.9245,1.7703,1,0.8266,0.5924,0,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Five extrema")
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')
hold on
%%
subplot(5,2,8)
x=0:0.01:10;
% alpha=1;
% s=1;
% k=1;
% nu1=2;
% nu2=1;
% nu3=-10;
params=[-0.3996,1.3606,1,0.7170,0.7436,0.7941,0];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=($-$0.3996,1.3606,1,0.7170,0.7436,0.7941,0)",LineWidth=2,color=[0.35 0.35 0.35])
title("Seven extrema")
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')
hold on
%%
subplot(5,2,9)
x=0:0.01:10;
% alpha=1;
% s=1;
% k=1;
% nu1=2;
% nu2=-1;
% nu3=-2;
% nu4=-2;
params=[-0.1704,1.1232,1,0.7111,0.6537,0.7122,0.6885];
[f,F,h,moments]=new_gled_prop(params,x);

plot(x,h, 'DisplayName',"\boldmath$\theta$=($-$0.1704,1.1232,1,0.7111,0.6537,0.7122,0.6885)",LineWidth=2,color=[0.35 0.35 0.35])
title("Nine extrema")
hl = legend('show','Location','northeast');
set(hl, 'Interpreter','latex')
saveas(gcf,"hazard.eps")