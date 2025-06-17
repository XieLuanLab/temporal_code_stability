%% OOO Components *drift* vs rate *drift* 
nVS=1
addpath(genpath(pwd));
stimNames = {'DG','SG','RFG','','IJO'};
nSamplesToTake=30;

ZeroMean=0;
DivideSum=1;
kernalWidth=24;
kfold=5;
perDayCVcomp=true;
usePCA=false;

load([pwd filesep stimNames{nVS} filesep 'VS_Store2025_selectedOnlyRefit.mat'])
load('nonNoise_whichAnimalLayerShankID7.mat')
load([pwd filesep stimNames{nVS} filesep 'pFinalFinalSelect.mat'])

cd([pwd filesep stimNames{nVS}])
VS_Store(avail==0)={nan(size(VS_Store{1}))};
VS_Store=cellfun(@(x) x',VS_Store,'UniformOutput',0);

    for ani=1:5
    filename=['LDA_' num2str(kernalWidth*2) 'ms_' num2str(kfold) 'pD' num2str(perDayCVcomp) 'CV_a_' num2str(ani) 'nVS_' num2str(nVS) '_ZeroMean_' num2str(ZeroMean) '_DivideSum_' num2str(DivideSum) '_PCA' num2str(usePCA) '.mat'];
    load(filename,['FR_Store'])
    VS_Store(whichAnimalPerUnit==ani,:)=FR_Store;
    VS_Store=fillEmptyAsNaN(VS_Store);
    end    

VS_Store=cellfun(@(x) x',VS_Store,'UniformOutput',0);

load('umapWorstDistance.mat');
umapThres=85;
umapWorstDistance=cat(1,maxInterDayDistStore{:});WAPU=whichAnimalPerUnit(finalSelect,:);
umapWorstDistanceCri=false(size(umapWorstDistance));
for ani=1:5
umapWorstDistanceCri(WAPU==ani)=umapWorstDistance(WAPU==ani)<prctile(umapWorstDistance(WAPU==ani),umapThres);
end

VS_Store=VS_Store(finalSelect,:);
testPara=VS_Store;VS_Store=[];



switch nVS
    case 1
        trialN=50
        stimN=16
    case 5
         trialN=60
        stimN=100
    case 2
        trialN=90
        stimN=30
        case 3
        trialN=30
        stimN=81
end
   
    testPara(~pFinal)={[]};
    testPara=fillEmptyAsNaN(testPara);
    testPara=testPara(umapWorstDistanceCri,:);

rng(1234)
nPercent=50;nSampledTrials=floor(trialN*nPercent/100);
takeNpercentSamples = @(x) mean(x(randperm(trialN,nSampledTrials),:));
takeOddSamples = @(x) mean(x(1:2:end,:));
takeEvenSamples = @(x) mean(x(2:2:end,:));

take1_3Samples = @(x) mean(x(1:3:end,:));
take2_3Samples = @(x) mean(x(2:3:end,:));
take3_3Samples = @(x) mean(x(3:3:end,:));


  driftRate = nan(size(testPara,1),size(testPara,2),nSamplesToTake);
  driftRate2 = nan(size(testPara,1),size(testPara,2));
  driftRate3 = nan(size(testPara,1),size(testPara,2));

  for i = 1:size(testPara,1)
    temp=testPara(i,:);

     tempOddSample=cellfun(@(x) takeOddSamples(x),temp,'UniformOutput',0);
    tempEvenSample=cellfun(@(x) takeEvenSamples(x),temp,'UniformOutput',0);
    tempSampleOddEven=[cat(1,tempOddSample{:}) cat(1,tempEvenSample{:})];
    tempSampleOddEvenDay2Day=corr(tempSampleOddEven');


 temp1_3Sample=cellfun(@(x) take1_3Samples(x),temp,'UniformOutput',0);
 temp2_3Sample=cellfun(@(x) take2_3Samples(x),temp,'UniformOutput',0);
  temp3_3Sample=cellfun(@(x) take3_3Samples(x),temp,'UniformOutput',0);
temp3=[cat(1,temp1_3Sample{:}) cat(1,temp2_3Sample{:}) cat(1,temp3_3Sample{:})];
temp3=corr(temp3');


     for d = 1:14
   driftRate2(i,d+1)=nanmean(diag(tempSampleOddEvenDay2Day,d));
    driftRate3(i,d+1)=nanmean(diag(temp3,d));
     end

    for s=1:nSamplesToTake
    tempSample=cellfun(@(x) takeNpercentSamples(x),temp,'UniformOutput',0);
    tempSample=cat(1,tempSample{:});
    tempSampleDay2Day=corr(tempSample');

    for d = 1:14
    driftRate(i,d+1,s)=nanmean(diag(tempSampleDay2Day,d));
    end

    [s1,s2]=takeNpercentSamplesAndTheOtherHalf(temp,trialN,nSampledTrials);
   s2=cat(1,s2{:});
   s1=cat(1,s1{:});
   driftRate(i,1,s)=nanmean(arrayfun(@(x) corr(s1(x,:)',s2(x,:)'),1:15));
  
    
    end


    end

nthComp=1:3;
score15=[];
axis15=[];

for ani = 1:5
    ani
filename=['LDA_' num2str(kernalWidth*2) 'ms_' num2str(kfold) 'pD' num2str(perDayCVcomp) 'CV_a_' num2str(ani) 'nVS_' num2str(nVS) '_ZeroMean_' num2str(ZeroMean) '_DivideSum_' num2str(DivideSum) '_PCA' num2str(usePCA) '.mat'];
load(filename,'LDA_15score_Store')

LDA_15score_Store=fillEmptyAsNaN(LDA_15score_Store);
LDA_Comp=cellfun(@(x) x(:,:,nthComp),LDA_15score_Store,'UniformOutput',false);
score15=[score15;LDA_Comp];

end

score15=score15(finalSelect,:);


score15(pFinal==0)={[]};
score15=fillEmptyAsNaN(score15);
score15=score15(umapWorstDistanceCri,:);


driftRateCP = nan(size(score15,1),size(score15,2),nSamplesToTake);
driftRateCP1 = nan(size(score15,1),size(score15,2),nSamplesToTake);
driftRateCP2 = nan(size(score15,1),size(score15,2),nSamplesToTake);
driftRateCP3 = nan(size(score15,1),size(score15,2),nSamplesToTake);


    for i = 1:size(score15,1)
    temp=cat(1,score15(i,:));

    if numel(nthComp)>=1
    temp=cellfun(@(x) permute(x,[2 1 3]),temp,'Uni',0);
    temp1=cellfun(@(x) x(:,:,1),temp,'UniformOutput',0);
    temp2=cellfun(@(x) x(:,:,2),temp,'UniformOutput',0);
    temp3=cellfun(@(x) x(:,:,3),temp,'UniformOutput',0);
    temp=cellfun(@(x) x(:,:),temp,'UniformOutput',0);
    end

    for s=1:nSamplesToTake
    tempSample=cellfun(@(x) takeNpercentSamples(x),temp,'UniformOutput',0);
    tempSample=cat(1,tempSample{:});
    tempSampleDay2Day=corr(tempSample');
    for d = 1:14
    driftRateCP(i,d+1,s)=nanmean(diag(tempSampleDay2Day,d));
    end

    [s1,s2]=takeNpercentSamplesAndTheOtherHalf(temp,trialN,nSampledTrials);
   s2=cat(1,s2{:});
   s1=cat(1,s1{:});
   driftRateCP(i,1,s)=nanmean(arrayfun(@(x) corr(s1(x,:)',s2(x,:)'),1:15));

    
    end

    for s=1:nSamplesToTake
    tempSample=cellfun(@(x) takeNpercentSamples(x),temp1,'UniformOutput',0);
    tempSample=cat(1,tempSample{:});
    tempSampleDay2Day=corr(tempSample');
    for d = 1:14
    driftRateCP1(i,d+1,s)=nanmean(diag(tempSampleDay2Day,d));
    end

    [s1,s2]=takeNpercentSamplesAndTheOtherHalf(temp1,trialN,nSampledTrials);
   s2=cat(1,s2{:});
   s1=cat(1,s1{:});
   driftRateCP1(i,1,s)=nanmean(arrayfun(@(x) corr(s1(x,:)',s2(x,:)'),1:15));

    end

    for s=1:nSamplesToTake
    tempSample=cellfun(@(x) takeNpercentSamples(x),temp2,'UniformOutput',0);
    tempSample=cat(1,tempSample{:});
    tempSampleDay2Day=corr(tempSample');
    for d = 1:14
    driftRateCP2(i,d+1,s)=nanmean(diag(tempSampleDay2Day,d));
    end

   [s1,s2]=takeNpercentSamplesAndTheOtherHalf(temp2,trialN,nSampledTrials);
   s2=cat(1,s2{:});
   s1=cat(1,s1{:});
   driftRateCP2(i,1,s)=nanmean(arrayfun(@(x) corr(s1(x,:)',s2(x,:)'),1:15));
   
    end

    for s=1:nSamplesToTake
    tempSample=cellfun(@(x) takeNpercentSamples(x),temp3,'UniformOutput',0);
    tempSample=cat(1,tempSample{:});
    tempSampleDay2Day=corr(tempSample');
    for d = 1:14
    driftRateCP3(i,d+1,s)=nanmean(diag(tempSampleDay2Day,d));
    end

   [s1,s2]=takeNpercentSamplesAndTheOtherHalf(temp3,trialN,nSampledTrials);
   s2=cat(1,s2{:});
   s1=cat(1,s1{:});
   driftRateCP3(i,1,s)=nanmean(arrayfun(@(x) corr(s1(x,:)',s2(x,:)'),1:15));

    end
    end
%%
driftRate(driftRate<0)=0;
driftRateCP(driftRateCP<0)=0;

dR=nanmean(driftRate,3);
dRCP=nanmean(driftRateCP,3);

driftRateMG=(dR(:,1)-dR(:,1:end))./(dR(:,1)+dR(:,1:end));
driftRateCPMG=(dRCP(:,1)-dRCP(:,1:end))./(dRCP(:,1)+dRCP(:,1:end));

h=figure('Position',[1144         221         300         405])
global xaxisshaded
xaxisshaded=0:14
hanlinErrorBar(driftRateMG,'b')
    hold on;
     hanlinErrorBar(driftRateCPMG,'r')
    xlabel('interval (days)')
    ylabel('representational drift index')
    xaxisshaded=[]
xlim([-0.5 14.5])
l=legend('rate','components','Color','none','EdgeColor','none','Location',[0.2616 0.8237 0.5222 0.1219])
set(gca,'FontSize',16)
l.FontSize=15
title([stimNames{nVS} '-all'],'FontWeight','normal')
if nVS==5
    title('NI-all','FontWeight','normal')
end
yl=ylim;

switch nVS
    case 1
ylim([0 0.176])
   case 2
ylim([0 0.145])
case 3
ylim([0 0.148])
 case 5
ylim([0 0.2070])
end
yl=ylim;

var1=driftRateCPMG
var2=driftRateMG
[r,c]=find(isnan(var1)~=isnan(var2))
valid=any(~isnan(var1),2) | any(~isnan(var2),2) ;

var1=var1(valid,:)
var2=var2(valid,:)
[nanmean(var1(:,end)) sem(var1(:,end))]
[nanmean(var2(:,end)) sem(var2(:,end))]

lme=hanlinLME([{var1} {var2}],0);
str1=sprintf('p(same slope)\n=%.2g',double(lme{4}.Coefficients(4,6)))

if double(lme{4}.Coefficients(4,6))<realmin
str1=sprintf('p(same slope)\n<%.2g',realmin)
end

text(0.25,yl(2)*0.81,str1,'FontSize',15,'Color','k')

