
if updatedRate==1
fr15=cell(5,1);
for ani = 1:5
filename=['LDA_' num2str(kernalWidth*2) 'ms_' num2str(kfold) 'pD' num2str(1) 'CV_a_' num2str(ani) 'nVS_' num2str(nVS) '_ZeroMean_' num2str(ZeroMean) '_DivideSum_' num2str(DivideSum) '_PCA' num2str(0) '.mat'];
load([basePath filesep stimNames{nVS} filesep filename],'FR_Store');
FR_Store=fillEmptyAsNaN(FR_Store);
fr15{ani,1}=FR_Store; 
end
VS_Store = cat(1,fr15{:,1});
VS_Store = cellfun(@(x) x',VS_Store,'UniformOutput',0);
end
load([basePath filesep stimNames{nVS} filesep 'pFinalFinalSelect.mat'])
load([basePath filesep 'nonNoise_whichAnimalLayerShankID7.mat'])
stimNames = {'DG','SG','RFG','','IJO'};
load('umapWorstDistance.mat');
umapWorstDistance=cat(1,maxInterDayDistStore{:});WAPU=whichAnimalPerUnit(finalSelect,:);
umapWorstDistanceCri=false(size(umapWorstDistance));
for ani=1:5
umapWorstDistanceCri(WAPU==ani)=umapWorstDistance(WAPU==ani)<prctile(umapWorstDistance(WAPU==ani),umapThres);
end
VS_Store=VS_Store(finalSelect,:);
VS_Store(pFinal==0)={nan*VS_Store{1}};
finalSelectIdx=find(finalSelect);
finalSelectIdx_animal=cell(5,1);
WAPU_valid =(whichAnimalPerUnit(finalSelect));
designMatrix = cell(5,15); % 5 ani 15 d 3 CP
for ani = 1:5
VS_1A15D=VS_Store(sum(pFinal,2)>=greaterThanNd &  WAPU_valid==ani & umapWorstDistanceCri,:);
finalSelectIdx_animal{ani}=finalSelectIdx(find(sum(pFinal,2)>=greaterThanNd &  WAPU_valid==ani));
 nUbynC_15D = cell2mat(cellfun(@mean,VS_1A15D,'UniformOutput',false));
 zmean15d=nanmean(nUbynC_15D ,2)';
 zstd15d=nanstd(nUbynC_15D ,[],2)';
    for d = 1:15
   VS_1A1D = VS_1A15D(:,d);
    nUbynC1 = cell2mat(cellfun(@(x) mean(x(1:sampleEveryNTrial:end,:)),VS_1A1D,'UniformOutput',false));
    nUbynC2 = cell2mat(cellfun(@(x) mean(x(2:sampleEveryNTrial:end,:)),VS_1A1D,'UniformOutput',false));
    nCbynU=[nUbynC1 nUbynC2]';
    if zmode==0
    designMatrix{ani,d}=nCbynU;
    end
    end
end

