clc
clear all
load('volume_all_acd_results.mat')
names = {'BUD','ADBE','AZN','LMT','SHW','SCCO'};
acds = {'E-ACD','W-ACD','GG-ACD','B-ACD','Glam$_1$-ACD','Glam$_2$-ACD','Glam$_3$-ACD','Glam$_4$-ACD',...
    'E-ACD$_1$','W-ACD$_1$','GG-ACD$_1$','B-ACD$_1$','Glam$_1$-ACD$_1$','Glam$_2$-ACD$_1$','Glam$_3$-ACD$_1$','Glam$_4$-ACD$_1$',...
    'E-ACD$_2$','W-ACD$_2$','GG-ACD$_2$','B-ACD$_2$','Glam$_1$-ACD$_2$','Glam$_2$-ACD$_2$','Glam$_3$-ACD$_2$','Glam$_4$-ACD$_2$'};
for item = 1:length(names)
    figure('Renderer', 'painters', 'Position', [10 10 1200 600])
    for i = 1:24
        subplot(6,4,i)
        [Lower_bound,Upper_bound]=generate_confidence_interval(outzs{item}(:,i));
        
        histogram(outzs{item}(:,i),20,'Normalization','probability','FaceColor','k')
        if item==1 && i==10
            ylim([0,1])
        elseif item==2 && i==10
            ylim([0,1])
        elseif item==3 && i==20
            ylim([0,1])
        elseif item==6 && i==10
            ylim([0,1])
        else
            ylim([0,0.11])
        end
        yline([Lower_bound,Upper_bound],'--')
        title(acds{i},'Interpreter','latex')
    end
    saveas(gcf,sprintf('./Figures/%s_outofsample_his',names{item}),'epsc')
end
%%
for item = 1:6
    figure('Renderer', 'painters', 'Position', [10 10 1200 600])
    for i = 1:24
        subplot(6,4,i)
        acf=autocorr(outzs{item}(:,i),NumLags=15);
        
        ylower=-1.96/sqrt(length(outzs{item}(:,i))-15);
        yupper=1.96/sqrt(length(outzs{item}(:,i))-15);
        
        yline([ylower,yupper],'-')
        hold on
        plot(acf(2:end),'o-','MarkerSize',3,color='k')
        xlabel('Lag')
        %ylim([floor(ylower*100)/100 ceil(acf(2)*10)/10])
        if item==2 && i==1
            ylim([-0.4,0.6])
        elseif item==2 && i==10
            ylim([-0.4,0.6])
        elseif item==2 && i==21
            ylim([-0.4,0.5])
        elseif item==2 && i==18
            ylim([-0.4,0.8])
        elseif item==3 && i==20
            ylim([-0.4,0.6])
        else
            ylim([-0.4,0.4])
        end
        title(acds{i},'Interpreter','latex')
    end
    saveas(gcf,sprintf('./Figures/%s_outofsample_autocorr',names{item}),'epsc')
end
