clc
clear all
load('volume_all_acd_results.mat')
names = {'BUD','ADBE','AZN','LMT','SHW','SCCO'};
figure('Renderer', 'painters', 'Position', [10 10 1200 600])
for item =1:6
    final_table = readtable(sprintf('./new2024/%s_volume_std.xlsx',names{item}));
    %% step 1
    sta_durations=final_table.std_duration;
    inx=sta_durations(1:round(length(sta_durations)*2/3-1),1);
    if item~=6
        parameters = params{item,9}(:,1);
        alpha = parameters(1);
        beta = parameters(2);
        gamma = parameters(3);
        psi = zeros(1,length(inx));
        psi(1) = alpha/(1-beta-gamma);
        
        for i=2:length(inx)
            psi(i) = alpha + beta*log(inx(i-1)) + gamma*psi(i-1);
        end
        epsilon = inx./exp(psi');
        if item==1
            parameters1 = params{item,15}(:,1);
            h1 = Glam_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
            subplot(3,2,item)
        % plot(t,h,'o','MarkerSize',3,color='k')
        % hold on
            plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
            ylabel("Hazard")
            legend({'Glam_3-ACD_1'},'Location','southeast')
            title(sprintf("%s",names{item}))
            hold on
        elseif item==2
            parameters1 = params{item,10}(:,1);
            h1 = weibull_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
            subplot(3,2,item)
            plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
            ylabel("Hazard")
            legend({'W-ACD_1'},'Location','southeast')
            title(sprintf("%s",names{item}))
            hold on
        elseif item==3
            parameters1 = params{item,12}(:,1);
            h1 = burr_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
            subplot(3,2,item)
            plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
            ylabel("Hazard")
            legend({'B-ACD_1'},'Location','southeast')
            title(sprintf("%s",names{item}))
            hold on
        elseif item==4
            parameters1 = params{item,14}(:,1);
            h1 = Glam_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
            subplot(3,2,item)
            plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
            ylabel("Hazard")
            legend({'Glam_2-ACD_1'},'Location','southeast')
            title(sprintf("%s",names{item}))
            hold on
        else
            parameters1 = params{item,16}(:,1);
            h1 = Glam_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
            subplot(3,2,item)
            plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
            ylabel("Hazard")
            legend({'Glam_4-ACD_1'},'Location','southeast')
            title(sprintf("%s",names{item}))
            hold on
        end
    else
        parameters = params{item,1}(:,1);
        alpha = parameters(1);
        beta = parameters(2);
        gamma = parameters(3);
        psi = zeros(1,length(inx));
        psi(1) = alpha/(1-beta-gamma);
        for i=2:length(inx)
            psi(i) = alpha + beta*inx(i-1) + gamma*psi(i-1);
        end
        epsilon=inx./psi';
        k=sum(inx>quantile(inx,0.95));
        [t,h] = baseline_hazard(epsilon,k);
        
        parameters1 = params{item,5}(:,1);
        h1 = Glam_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
        subplot(3,2,item)
        % plot(t,h,'o','MarkerSize',3,color='k')
        % hold on
        plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
        ylabel("Hazard")
        legend({'Glam_1-ACD'},'Location','southeast')
        title(sprintf("%s",names{item}))
    end
    %     % hold on
    %     % parameters = params{item,2}(:,1);
    %     % clear gamma
    %     % plot(0.01:0.01:max(t),parameters(end)*gamma(1+1/parameters(end))*((0.01:0.01:max(t)+5)*gamma(1+1/parameters(end))).^(parameters(end)-1))
    % else
    %     parameters = params{item,9}(:,1);
    %     alpha = parameters(1);
    %     beta = parameters(2);
    %     gamma = parameters(3);
    %     psi = zeros(1,length(inx));
    %     psi(1) = alpha/(1-beta-gamma);
    % 
    %     for i=2:length(inx)
    %         psi(i) = alpha + beta*log(inx(i-1)) + gamma*psi(i-1);
    %     end
    %     epsilon = inx./exp(psi');
    %     k=sum(inx>quantile(inx,0.95));
    %     [t,h] = baseline_hazard(epsilon,k);
    %     parameters1 = params{item,16}(:,1);
    %     h1 = Glam_hazard(min(epsilon):0.01:max(epsilon),parameters1(4:end)');
    %     subplot(3,2,item)
    %     % plot(t,h,'o','MarkerSize',3,color='k')
    %     % hold on
    %     plot(min(epsilon):0.01:max(epsilon),h1,LineWidth=1.5,Color='k')
    %     ylabel("Hazard")
    %     legend({'Glam_4-ACD_1'},'Location','best')
    %     title(sprintf("%s",names{item}))
    %     hold on
    % 
    % end
end