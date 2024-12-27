%%
clc
clear all
load('volume_all_acd_results.mat')

llfs = final_llf([1,2,3,4,5,6],:)'

aics = [-2;-2;-2;-2;2;2;2;2;-2;-2;-2;-2;2;2;2;2;-2;-2;-2;-2;2;2;2;2].*llfs+2*[3;4;5;5;4;5;6;7;3;4;5;5;4;5;6;7;3;4;5;5;4;5;6;7]
bics = [-2;-2;-2;-2;2;2;2;2;-2;-2;-2;-2;2;2;2;2;-2;-2;-2;-2;2;2;2;2].*llfs+[3;4;5;5;4;5;6;7;3;4;5;5;4;5;6;7;3;4;5;5;4;5;6;7]*[log(843),...
    log(241),log(463),log(371),log(707),log(638)]
%%
clc
clear all
load('volume_all_acd_results.mat')
for item=1:6
    item
    [(1:24)',test_results{1,item}(:,[5,1,2])]
    [(1:24)',test_results{1,item}(:,[6,3,4])]
end
%%
clc
clear all
load('volume_linear_acd_results_data8.mat')
llfs = final_llf([3,14,13,20,27,30],:)'
aics = [-2;-2;-2;-2;2;2;2;2].*final_llf([3,14,13,20,27,30],:)'+2*[3;4;5;5;4;5;6;7]
bics = [-2;-2;-2;-2;2;2;2;2].*final_llf([3,14,13,20,27,30],:)'+[3;4;5;5;4;5;6;7]*[log(843),...
    log(241),log(463),log(371),log(707),log(638)]


%%
clc
clear all
load('volume_log_acd1_results_data8.mat')
llfs = final_llf([3,14,13,20,27,30],:)'
aics = [-2;-2;-2;-2;2;2;2;2].*final_llf([3,14,13,20,27,30],:)'+2*[3;4;5;5;4;5;6;7]
bics = [-2;-2;-2;-2;2;2;2;2].*final_llf([3,14,13,20,27,30],:)'+[3;4;5;5;4;5;6;7]*[log(843),...
    log(241),log(463),log(371),log(707),log(638)]
%%
clc
clear all
load('volume_log_acd2_results_data8.mat')
llfs = final_llf([3,14,13,20,27,30],:)'
aics = [-2;-2;-2;-2;2;2;2;2].*final_llf([3,14,13,20,27,30],:)'+2*[3;4;5;5;4;5;6;7]
bics = [-2;-2;-2;-2;2;2;2;2].*final_llf([3,14,13,20,27,30],:)'+[3;4;5;5;4;5;6;7]*[log(843),...
    log(241),log(463),log(371),log(707),log(638)]
%% out of sample 
clc
clear all
load('volume_linear_acd_results_data8.mat')
results = [];
for i=13%[7,9,13,15,17,25]%[3,14,24,20,27,30]
    results = [results,[test_results{1,i}(:,end),test_results{1,i}(:,3:4)]];
end
results
%% out of sample 
clc
clear all
load('volume_log_acd1_results_data8.mat')
results = [];
for i=13%[7,9,13,15,17,25]%[3,14,24,20,27,30]
    results = [results,[test_results{1,i}(:,end),test_results{1,i}(:,3:4)]];
end
results
%% out of sample 
clc
clear all
load('volume_log_acd2_results_data8.mat')
results = [];
for i=13%[7,9,13,15,17,25]%[3,14,24,20,27,30]
    results = [results,[test_results{1,i}(:,end),test_results{1,i}(:,3:4)]];
end
results


%% in sample 
clc
clear all
load('volume_linear_acd_results_data8.mat')
results = [];
for i=[3,14,13,20,27,30]%9%[7,9,13,15,17,25]%[3,14,24,20,27,30]
    results = [results,[test_results{1,i}(:,end-1),test_results{1,i}(:,1:2)]];
end
results
%% in sample 
clc
clear all
load('volume_log_acd1_results_data8.mat')
results = [];
for i=[3,14,13,20,27,30]%9%[7,9,13,15,17,25]%[3,14,24,20,27,30]
    results = [results,[test_results{1,i}(:,end-1),test_results{1,i}(:,1:2)]];
end
results
%% in sample 
clc
clear all
load('volume_log_acd2_results_data8.mat')
results = [];
for i=[3,14,13,20,27,30]%9%[7,9,13,15,17,25]%[3,14,24,20,27,30]
    results = [results,[test_results{1,i}(:,end-1),test_results{1,i}(:,1:2)]];
end
results
