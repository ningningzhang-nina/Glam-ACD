clc
clear all
names = {'BUD','ADBE','AZN','LMT','SHW','SCCO'};
%names = {'CSCO','RY','AZN','AMAT','SBUX','HCA'};
%names={'RY','SBUX'};
lines = {'-o','--v','-d','-s','-p','-h'};
figure('Renderer', 'painters', 'Position', [10 10 1200 600])
for item = 1:length(names)
    final_table = readtable(sprintf('./new2024/%s_volume_std.xlsx',names{item}));
    x = final_table.std_duration;
    s = sqrt(var(x));
    n = length(x);
    b = (0.9*s*n)^(-0.4);
    pts = min(x):0.05:max(x);
    ft = gamma_ksdensity(x,pts,b);
    % figure
    plot(pts,ft,lines{item},"LineWidth",2.5,color=[0.35 0.35 0.35])
    %title(sprintf('%s',names{item}))
    hold on
    [h,pValue,stat,cValue]=lbqtest(final_table.duration,Lags=10);
    names{item}
    [h,pValue]
    [length(x),min(x),max(x),mean(x),sqrt(var(final_table.duration))/mean(final_table.duration),stat]
    min(x)
    % figure
    % ksdensity(x)
    % title(sprintf('%s_volume_density',names{item}))
    
    % [acf,lags]=autocorr(final_table.std_duration,NumLags=50);
    % plot(1:50,acf(2:end),lines{item},"LineWidth",2.5)
    % ylabel('Autocorrelation')
    % xlabel('Lags')
    % hold on
    %title(sprintf('%s',names{item}))
    %saveas(gcf,sprintf('./Figures/%s_acf.eps',names{item}))
    % figure
    % autocorr(x)
    
end
legend({'BUD','ADBE','AZN','LMT','SHW','SCCO'})
%legend({'CSCO','RY','AZN','AMAT','SBUX','HCA'})
saveas(gcf,'./Figures/density_volume_durations.eps','epsc')
%%
figure('Renderer', 'painters', 'Position', [10 10 1200 600])
for item = 1:length(names)
    final_table = readtable(sprintf('./new2024/%s_volume_std.xlsx',names{item}));
    x = final_table.std_duration;
    s = sqrt(var(x));
    n = length(x);
    b = (0.9*s*n)^(-0.4);
    pts = min(x):0.05:max(x);
    % ft = gamma_ksdensity(x,pts,b);
    % % figure
    % plot(pts,ft,lines{item},"LineWidth",2.5,color=[0.35 0.35 0.35])
    % %title(sprintf('%s',names{item}))
    % hold on
    % [h,pValue,stat,cValue]=lbqtest(final_table.duration,Lags=10);
    % names{item}
    % [h,pValue]
    % [length(x),min(x),max(x),mean(x),sqrt(var(final_table.duration))/mean(final_table.duration),stat]
    % min(x)
    % figure
    % ksdensity(x)
    % title(sprintf('%s_volume_density',names{item}))
    
    [acf,lags]=autocorr(final_table.std_duration,NumLags=50);
    plot(1:50,acf(2:end),lines{item},"LineWidth",2.5,color=[0.35 0.35 0.35])
    ylabel('Autocorrelation')
    xlabel('Lags')
    hold on
    %title(sprintf('%s',names{item}))
    %saveas(gcf,sprintf('./Figures/%s_acf.eps',names{item}))
    % figure
    % autocorr(x)
    
end
legend({'CSCO','RY','AZN','AMAT','SBUX','HCA'})
%legend({'BUD','ADBE','CI','LMT','SHW','SCCO'})
%saveas(gcf,'./Figures/acf','epsc')
%saveas(gcf,'./Figures/density_volume_durations.eps','epsc')