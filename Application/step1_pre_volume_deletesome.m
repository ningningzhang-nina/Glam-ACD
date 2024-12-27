clear all
%names = {'KKR','TD','BMY','IBN','MDT','AMT','UBS','CI','BP','INFY','HCA','SO','SMFG','SHW','ICE','DUK','SCCO','COP','BLK','TJX','SYK','PGR','LMT'};
names = {'MRK','CVX','TM','CRM','SAP','TMO','NVS','ACN','DHR','GE','ABT'};
for item =length(names)
T = readtable(sprintf('/Users/gina/Downloads/%s.csv',names{item}));%30971766
%% step 1: 只保留9:30am to 4:00pm的数据
T.EX=char(T.EX);
T=T((T.TIME_M>='09:30:00')&(T.TIME_M<='16:00:00')&((T.EX=='N')|(T.EX=='A')|(T.EX=='P')|(T.EX=='C')|(T.EX=='M')),:); % 392962 条数据
%% step 2: 同一个时间戳进行合并
TT=groupsummary(T,["DATE",'TIME_M'],["sum"],"SIZE"); %  392962 条数据
%%
days = unique(TT.DATE);
TTT = TT(TT.DATE==days(1),:);
if TTT(2,:).sum_SIZE>1000
    TTT(2,:) = [];
end
if TTT(1,:).sum_SIZE>1000
    TTT(1,:) = [];
end
for day =2:length(days)
    a=TT(TT.DATE==days(day),:);
    if a(2,:).sum_SIZE>1000
        a(2,:) = [];
    end
    if a(1,:).sum_SIZE>1000
        a(1,:) = [];
    end
    TTT = [TTT;a];
end
%% step 2: open quote data
T1 = readtable(sprintf('/Users/gina/Downloads/%s1.csv',names{item})); 
% step 1: 只保留9:30am to 4:00pm的数据
T1.EX=char(T1.EX);
T1=T1((T1.TIME_M>='09:30:00')&(T1.TIME_M<='16:00:00')&((T1.EX=='N')|(T1.EX=='A')|(T1.EX=='P')|(T1.EX=='C')|(T1.EX=='M')),:); % 23789125 条数据

%%
%TTT(TTT.DATE==days(80),:)=[];
days = unique(TTT.DATE);
%% step 3: 初步得到volume duration,这里考虑的两个时间点之间的volume是只会把前面的volume计算进去，最后一个时间点的volume不计算进去
criterion = 25000;
for i = 1:length(days)
    TTTT = TTT(TTT.DATE==days(i),:);
    temp_T = T1(T1.DATE==days(i),:);
    volumetable = temp_T(1,{'DATE','TIME_M'});
    j = 1;
    n = size(TTTT,1);
    while j<n
        cumvolume=cumsum(TTTT(j:end,:).sum_SIZE);
        locs = find(cumvolume>=criterion);
        if isempty(locs)
            break
        else
            j = locs(1)+j;
            if j<=n
                try
                    find_next_time = temp_T(temp_T.TIME_M>TTTT.TIME_M(j),:);
                    volumetable=[volumetable;find_next_time(1,{'DATE','TIME_M'})]; 
                catch
                    break
                end
            else
                break
            end
        end   
    end
    volume_duration = diff(volumetable.TIME_M);
    [h,m,s]=hms(volume_duration);
    volume_duration = h*60*60+60*m+s;
    if i==1
        volumetableall = volumetable(1:end-1,:);
        volumetableall.duration = volume_duration;
    else
        volumetable.duration = [volume_duration;0];
        volumetableall = [volumetableall;volumetable(1:end-1,:)];
 
    end
end
[h,m,s]=hms(volumetableall.TIME_M);
volumetableall.H=h;
volumetableall.M=m;
volumetableall.S=s;

%% step 4: 去除volumetableallily 的趋势用cubic spline
knots = {'09:30:00';'10:00:00';'10:30:00';'11:00:00';'11:30:00';'12:00:00';'12:30:00';'13:00:00';...
    '13:30:00';'14:00:00';'14:30:00';'15:00:00';'15:30:00';'16:00:00'};
volumetableall.ind1 = (volumetableall.TIME_M<knots(2))&(volumetableall.TIME_M>=knots(1));
volumetableall.ind2 = (volumetableall.TIME_M<knots(3))&(volumetableall.TIME_M>=knots(2));
volumetableall.ind3 = (volumetableall.TIME_M<knots(4))&(volumetableall.TIME_M>=knots(3));
volumetableall.ind4 = (volumetableall.TIME_M<knots(5))&(volumetableall.TIME_M>=knots(4));
volumetableall.ind5 = (volumetableall.TIME_M<knots(6))&(volumetableall.TIME_M>=knots(5));
volumetableall.ind6 = (volumetableall.TIME_M<knots(7))&(volumetableall.TIME_M>=knots(6));
volumetableall.ind7 = (volumetableall.TIME_M<knots(8))&(volumetableall.TIME_M>=knots(7));
volumetableall.ind8 = (volumetableall.TIME_M<knots(9))&(volumetableall.TIME_M>=knots(8));
volumetableall.ind9 = (volumetableall.TIME_M<knots(10))&(volumetableall.TIME_M>=knots(9));
volumetableall.ind10 = (volumetableall.TIME_M<knots(11))&(volumetableall.TIME_M>=knots(10));
volumetableall.ind11 = (volumetableall.TIME_M<knots(12))&(volumetableall.TIME_M>=knots(11));
volumetableall.ind12 = (volumetableall.TIME_M<knots(13))&(volumetableall.TIME_M>=knots(12));
volumetableall.ind13 = (volumetableall.TIME_M<knots(14))&(volumetableall.TIME_M>=knots(13));

volumetableall.knot_duration = ones(size(volumetableall,1),1);
for i=1:13
    [~,m,s]=hms(volumetableall(eval(strcat('volumetableall.ind', string(i))) == 1,:).TIME_M-duration(knots(i)));
    volumetableall(eval(strcat('volumetableall.ind', string(i))) == 1,:).knot_duration=m*60+s;
end

%% least square with constraints
XX = [volumetableall.ind1,volumetableall.ind1.*volumetableall.knot_duration,volumetableall.ind1.*volumetableall.knot_duration.^2,volumetableall.ind1.*volumetableall.knot_duration.^3,...
    volumetableall.ind2,volumetableall.ind2.*volumetableall.knot_duration,volumetableall.ind2.*volumetableall.knot_duration.^2,volumetableall.ind2.*volumetableall.knot_duration.^3,...
    volumetableall.ind3,volumetableall.ind3.*volumetableall.knot_duration,volumetableall.ind3.*volumetableall.knot_duration.^2,volumetableall.ind3.*volumetableall.knot_duration.^3,...
    volumetableall.ind4,volumetableall.ind4.*volumetableall.knot_duration,volumetableall.ind4.*volumetableall.knot_duration.^2,volumetableall.ind4.*volumetableall.knot_duration.^3,...
    volumetableall.ind5,volumetableall.ind5.*volumetableall.knot_duration,volumetableall.ind5.*volumetableall.knot_duration.^2,volumetableall.ind5.*volumetableall.knot_duration.^3,...
    volumetableall.ind6,volumetableall.ind6.*volumetableall.knot_duration,volumetableall.ind6.*volumetableall.knot_duration.^2,volumetableall.ind6.*volumetableall.knot_duration.^3,...
    volumetableall.ind7,volumetableall.ind7.*volumetableall.knot_duration,volumetableall.ind7.*volumetableall.knot_duration.^2,volumetableall.ind7.*volumetableall.knot_duration.^3,...
    volumetableall.ind8,volumetableall.ind8.*volumetableall.knot_duration,volumetableall.ind8.*volumetableall.knot_duration.^2,volumetableall.ind8.*volumetableall.knot_duration.^3,...
    volumetableall.ind9,volumetableall.ind9.*volumetableall.knot_duration,volumetableall.ind9.*volumetableall.knot_duration.^2,volumetableall.ind9.*volumetableall.knot_duration.^3,...
    volumetableall.ind10,volumetableall.ind10.*volumetableall.knot_duration,volumetableall.ind10.*volumetableall.knot_duration.^2,volumetableall.ind10.*volumetableall.knot_duration.^3,...
    volumetableall.ind11,volumetableall.ind11.*volumetableall.knot_duration,volumetableall.ind11.*volumetableall.knot_duration.^2,volumetableall.ind11.*volumetableall.knot_duration.^3,...
    volumetableall.ind12,volumetableall.ind12.*volumetableall.knot_duration,volumetableall.ind12.*volumetableall.knot_duration.^2,volumetableall.ind12.*volumetableall.knot_duration.^3,...
    volumetableall.ind13,volumetableall.ind13.*volumetableall.knot_duration,volumetableall.ind13.*volumetableall.knot_duration.^2,volumetableall.ind13.*volumetableall.knot_duration.^3];
Aeq = [1 zeros(1,4*13-1);...
    1 180 180^2 180^3 -1 zeros(1,4*12-1);...
    zeros(1,4) 1 180 180^2 180^3 -1 zeros(1,4*11-1);...
    zeros(1,4*2) 1 180 180^2 180^3 -1 zeros(1,4*10-1);...
    zeros(1,4*3) 1 180 180^2 180^3 -1 zeros(1,4*9-1);...
    zeros(1,4*4) 1 180 180^2 180^3 -1 zeros(1,4*8-1);...
    zeros(1,4*5) 1 180 180^2 180^3 -1 zeros(1,4*7-1);...
    zeros(1,4*6) 1 180 180^2 180^3 -1 zeros(1,4*6-1);...
    zeros(1,4*7) 1 180 180^2 180^3 -1 zeros(1,4*5-1);...
    zeros(1,4*8) 1 180 180^2 180^3 -1 zeros(1,4*4-1);...
    zeros(1,4*9) 1 180 180^2 180^3 -1 zeros(1,4*3-1);...
    zeros(1,4*10) 1 180 180^2 180^3 -1 zeros(1,4*2-1);...
    zeros(1,4*11) 1 180 180^2 180^3 -1 zeros(1,4*1-1);...
    0 1 2*180 3*180^2 0 -1 zeros(1,4*12-2);...
    zeros(1,4) 0 1 2*180 3*180^2 0 -1 zeros(1,4*11-2);...
    zeros(1,4*2) 0 1 2*180 3*180^2 0 -1 zeros(1,4*10-2);...
    zeros(1,4*3) 0 1 2*180 3*180^2 0 -1 zeros(1,4*9-2);...
    zeros(1,4*4) 0 1 2*180 3*180^2 0 -1 zeros(1,4*8-2);...
    zeros(1,4*5) 0 1 2*180 3*180^2 0 -1 zeros(1,4*7-2);...
    zeros(1,4*6) 0 1 2*180 3*180^2 0 -1 zeros(1,4*6-2);...
    zeros(1,4*7) 0 1 2*180 3*180^2 0 -1 zeros(1,4*5-2);...
    zeros(1,4*8) 0 1 2*180 3*180^2 0 -1 zeros(1,4*4-2);...
    zeros(1,4*9) 0 1 2*180 3*180^2 0 -1 zeros(1,4*3-2);...
    zeros(1,4*10) 0 1 2*180 3*180^2 0 -1 zeros(1,4*2-2);...
    zeros(1,4*11) 0 1 2*180 3*180^2 0 -1 zeros(1,4*1-2)];
beq = [mean(volumetableall.duration);zeros(24,1)];
x = lsqlin(XX,volumetableall.duration,[],[],Aeq,beq,[],[])
volumetableall.std_duration = volumetableall.duration./(XX*x);
%%
names{item}
min(volumetableall.std_duration)
%%
figure
ksdensity(volumetableall.std_duration)
xlim([0,6])
title('volume duration')
saveas(gcf,sprintf('./new2024/%s_std_volume2024_ksdensity.png',names{item}))
%% step 6: save data
writetable(volumetableall,sprintf('./new2024/%s_volume2024_std.xlsx',names{item}),'PreserveFormat',false);
end