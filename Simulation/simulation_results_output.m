clc
clear all
load('simulation_glgd1_acd_data1_n500_r2000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0356    0.0440   -0.0088   -0.0028
sd=std(params5(1:500,:))
% 0.0380    0.0426    0.0452    0.0538
ad=mean(sd5(1:500,:))
% 0.0377    0.0405    0.0445    0.0487
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_acd_data1_n1000_r2000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0269    0.0470   -0.0054   -0.0006
sd=std(params5(1:500,:))
% 0.0254    0.0297    0.0314    0.0354
ad=mean(sd5(1:500,:))
% 0.0248    0.0286    0.0306    0.0342
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd1_data1_n500_r1000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0462   -0.0026   -0.0098   -0.0069
sd=std(params5(1:500,:))
% 0.0168    0.0317    0.0779    0.0536
ad=mean(sd5(1:500,:))
% 0.0166    0.0291    0.0778    0.0483
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd1_data1_n1000_r1000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0457   -0.0007   -0.0058   -0.0032
sd=std(params5(1:500,:))
% 0.0119    0.0228    0.0603    0.0352
ad=mean(sd5(1:500,:))
% 0.0118    0.0210    0.0566    0.0340
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd2_data1_n500_r1000.mat')
mean(params5([65:564],:)-[alpha,beta,gamma,phi01],1)
% 0.0273    0.0409   -0.0083   -0.0096
sd=std(params5([65:564],:))
% 0.0404    0.0412    0.0361    0.0623
ad=mean(sd5([65:564],:))
% 0.0368    0.0388    0.0369    0.0481
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd2_data1_n1000_r1000.mat')
mean(params5([65:128,193:512,577:692],:)-[alpha,beta,gamma,phi01],1)
% 0.0239    0.0431   -0.0056   -0.0035
sd=std(params5([65:128,193:512,577:692],:))
% 0.0265    0.0277    0.0271    0.0389
ad=mean(sd5([65:128,193:512,577:692],:))
% 0.0258    0.0274    0.0258    0.0339
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_acd_data2_n500_r1000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0830    0.0770   -0.0081    0.0100
sd=std(params5(1:500,:))
% 0.0705    0.0391    0.0543    0.1775
ad=mean(sd5(1:500,:))
% 0.0703    0.0369    0.0532    0.1789
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_acd_data2_n1000_r1000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0743    0.0774   -0.0045    0.0062
sd=std(params5(1:500,:))
% 0.0486    0.0275    0.0389    0.1243
ad=mean(sd5(1:500,:))
% 0.0477    0.0259    0.0368    0.1264
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd1_data2_n500_r1000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0959   -0.0008   -0.0021    0.0087
sd=std(params5(1:500,:))
% 0.0391    0.0200    0.1088    0.1782
ad=mean(sd5(1:500,:))
% 0.0364    0.0182    0.1015    0.1795
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd1_data2_n1000_r1000.mat')
mean(params5(1:500,:)-[alpha,beta,gamma,phi01],1)
% 0.0957   -0.0004   -0.0025    0.0073
sd=std(params5(1:500,:))
% 0.0308    0.0141    0.0899    0.1254
ad=mean(sd5(1:500,:))
% 0.0297    0.0131    0.0842    0.1265
[alpha,beta,gamma,phi01]
%%
clc
clear all
load('simulation_glgd1_logacd2_data2_n500_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
mean(final_params5([1,3:16,18:23,25:66,68:127,129:286,289:291,293:389,391:406,408:414,416:437,439:468,470:500,502:506,508:513],:)-[alpha,beta,gamma,phi01],1)
% 0.0866    0.0667   -0.0019    0.0160
sd=std(final_params5([1,3:16,18:23,25:66,68:127,129:286,289:291,293:389,391:406,408:414,416:437,439:468,470:500,502:506,508:513],:))
% 0.0439    0.0402    0.0589    0.1969
ad=mean(final_sd5([1,3:16,18:23,25:66,68:127,129:286,289:291,293:389,391:406,408:414,416:437,439:468,470:500,502:506,508:513],:))
% 0.0422    0.0375    0.0613    0.1786
[alpha,beta,gamma,phi01]

%%
clc
clear all
load('simulation_glgd1_logacd2_data2_n1000_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
mean(final_params5([1:38,40:58,60:154,156:168,170:192,194:199,201:230,232:252,254:262,264:295,297:312,314:317,319:322,324:367,370:428,430:441,443:447,451:518,520:521],:)-[alpha,beta,gamma,phi01],1)
% 0.0831    0.0708   -0.0045    0.0079
sd=std(final_params5([1:38,40:58,60:154,156:168,170:192,194:199,201:230,232:252,254:262,264:295,297:312,314:317,319:322,324:367,370:428,430:441,443:447,451:518,520:521],:))
% 0.0324    0.0285    0.0455    0.1279
ad=mean(final_sd5([1:38,40:58,60:154,156:168,170:192,194:199,201:230,232:252,254:262,264:295,297:312,314:317,319:322,324:367,370:428,430:441,443:447,451:518,520:521],:))
% 0.0296    0.0264    0.0442    0.1263
[alpha,beta,gamma,phi01]

%%
clc
clear all
load('simulation_glgd2_acd_data1_n500_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
final_final_params5=final_params5((final_params5(:,4)<0)&(final_params5(:,5)<0),:);
final_final_sd5=final_sd5((final_params5(:,4)<0)&(final_params5(:,5)<0),:);
mean(final_final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0277    0.0380   -0.0068   -0.0146    0.0070
sd=std(final_final_params5(1:500,:))
% 0.0366    0.0420    0.0479    0.0765    0.0498
ad=mean(final_final_sd5(1:500,:))
% 0.0308    0.0373    0.0430    0.0660    0.0417
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_acd_data1_n1000_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
final_final_params5=final_params5((final_params5(:,4)<0)&(final_params5(:,5)<0),:);
final_final_sd5=final_sd5((final_params5(:,4)<0)&(final_params5(:,5)<0),:);
mean(final_final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0210    0.0403   -0.0040   -0.0103    0.0058
sd=std(final_final_params5(1:500,:))
% 0.0247    0.0289    0.0345    0.0550    0.0366
ad=mean(final_final_sd5(1:500,:))
% 0.0201    0.0267    0.0299    0.0477    0.0301
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd1_data1_n500_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
final_final_params5=final_params5((final_params5(:,end)<0)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
final_final_sd5=final_sd5((final_params5(:,end)<0)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
mean(final_final_params5-[alpha,beta,gamma,phi01,phi02],1)
% 0.0488    0.0120   -0.0380   -0.0186    0.0108
sd=std(final_final_params5)
% 0.0188    0.0249    0.0929    0.0795    0.0525
ad=mean(final_final_sd5)
% 0.0213    0.0220    0.1146    0.0713    0.0460
[alpha,beta,gamma,phi01,phi02]


%%
clc
clear all
load('simulation_glgd2_logacd1_data1_n1000_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
final_final_params5=final_params5((final_params5(:,end)<0)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
final_final_sd5=final_sd5((final_params5(:,end)<0)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
mean(final_final_params5-[alpha,beta,gamma,phi01,phi02],1)
% 0.0490    0.0073   -0.0372   -0.0064    0.0043
sd=std(final_final_params5)
% 0.0163    0.0185    0.0860    0.0555    0.0378
ad=mean(final_final_sd5)
% 0.0153    0.0144    0.0876    0.0493    0.0311
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd2_data1_n500_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
final_final_params5=final_params5((final_params5(:,end)<0)&(final_params5(:,end)>-2)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
final_final_sd5=final_sd5((final_params5(:,end)<0)&(final_params5(:,end)>-2)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
mean(final_final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0207    0.0392   -0.0166   -0.0035    0.0015
sd=std(final_final_params5(1:500,:))
% 0.0422    0.0436    0.0586    0.0849    0.0630
ad=mean(final_final_sd5(1:500,:))
% 0.0361    0.0383    0.0332    0.0755    0.0489
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd2_data1_n1000_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
final_final_params5=final_params5((final_params5(:,end)<-0.2)&(final_params5(:,end)>-2)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
final_final_sd5=final_sd5((final_params5(:,end)<-0.2)&(final_params5(:,end)>-2)&(final_params5(:,3)>0.5)&(final_params5(:,3)<1),:);
mean(final_final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0178    0.0396   -0.0106    0.0064   -0.0074
sd=std(final_final_params5(1:500,:))
% 0.0281    0.0303    0.0433    0.0542    0.0386
ad=mean(final_final_sd5(1:500,:))
% 0.0267    0.0285    0.0234    0.0538    0.0347
[alpha,beta,gamma,phi01,phi02]

%%
clc
clear all
load('simulation_glgd2_acd_data2_n500_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
mean(final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0362    0.0331   -0.0083   -0.0029    0.0018
sd=std(final_params5(1:500,:))
% 0.0551    0.0382    0.0585    0.0551    0.0680
ad=mean(final_sd5(1:500,:))
% 0.0468    0.0331    0.0504    0.0513    0.0624
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_acd_data2_n1000_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
mean(final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0313    0.0332   -0.0054   -0.0006    0.0016
sd=std(final_params5(1:500,:))
% 0.0381    0.0260    0.0415    0.0392    0.0471
ad=mean(final_sd5(1:500,:))
% 0.0300    0.0229    0.0336    0.0361    0.0438

[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd1_data2_n500_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
mean(final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0425    0.0048   -0.0164    0.0010    0.0022
sd=std(final_params5(1:500,:))
% 0.0143    0.0191    0.0368    0.0550    0.0665
ad=mean(final_sd5(1:500,:))
% 0.0112    0.0158    0.0315    0.0503    0.0591
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd1_data2_n1000_r1000.mat')
final_params5= [];
final_sd5 = [];
for i=1:1000
    if all(real(sd5(i,:)))
        final_params5 = [final_params5;params5(i,:)];
        final_sd5 = [final_sd5;sd5(i,:)];
    end
end
mean(final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0404    0.0026   -0.0083    0.0013   -0.0008
sd=std(final_params5(1:500,:))
% 0.0093    0.0129    0.0253    0.0382    0.0459
ad=mean(final_sd5(1:500,:))
% 0.0073    0.0107    0.0210    0.0357    0.0426
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd2_data2_n500_r1000.mat')

final_params5=params5((params5(:,end)<1)&(params5(:,end)>0)&(params5(:,3)>0.5)&(params5(:,4)<0)&(params5(:,4)>-1.5),:);
final_sd5=sd5((params5(:,end)<1)&(params5(:,end)>0)&(params5(:,3)>0.5)&(params5(:,4)<0)&(params5(:,4)>-1.5),:);
mean(final_params5(1:500,:)-[alpha,beta,gamma,phi01,phi02],1)
% 0.0203    0.0490   -0.0375    0.0238   -0.0018
sd=std(final_params5(1:500,:))
% 0.0368    0.0363    0.0641    0.0527    0.0651
ad=mean(final_sd5(1:500,:))
% 0.0368    0.0373    0.0709    0.0516    0.0627
[alpha,beta,gamma,phi01,phi02]
%%
clc
clear all
load('simulation_glgd2_logacd2_data2_n1000_r1000.mat')

final_params5=params5((params5(:,end)<1)&(params5(:,end)>0)&(params5(:,4)<0)&(params5(:,4)>-1.5),:);
final_sd5=sd5((params5(:,end)<1)&(params5(:,end)>0)&(params5(:,4)<0)&(params5(:,4)>-1.5),:);
mean(final_params5-[alpha,beta,gamma,phi01,phi02],1)
% 0.0274    0.0398   -0.0193    0.0133   -0.0015
sd=std(final_params5)
% 0.0252    0.0261    0.0482    0.0369    0.0420
ad=mean(final_sd5)
% 0.0253    0.0263    0.0489    0.0366    0.0443
[alpha,beta,gamma,phi01,phi02]