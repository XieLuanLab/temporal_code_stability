

load([basePath filesep stimNames{nVS}  filesep 'pFinalFinalSelect.mat'])
load('nonNoise_whichAnimalLayerShankID7.mat')
load('umapWorstDistance.mat');
umapWorstDistance=cat(1,maxInterDayDistStore{:});WAPU=whichAnimalPerUnit(finalSelect,:);
umapWorstDistanceCri=false(size(umapWorstDistance));
for ani=1:5
umapWorstDistanceCri(WAPU==ani)=umapWorstDistance(WAPU==ani)<prctile(umapWorstDistance(WAPU==ani),umapThres);
end

score15=[];
for ani = 1:5
    ani
   if  updatedComponents==1
filename=['LDA_' num2str(kernalWidth*2) 'ms_' num2str(kfold) 'pD' num2str(1) 'CV_a_' num2str(ani) 'nVS_' num2str(nVS) '_ZeroMean_' num2str(ZeroMean) '_DivideSum_' num2str(DivideSum) '_PCA' num2str(0) '.mat'];
   end
load([basePath filesep stimNames{nVS} filesep filename],'LDA_15score_Store');% Apr7 Submission version. 
LDA_15score_Store=fillEmptyAsNaN(LDA_15score_Store);
LDA_Comp=cellfun(@(x) x(:,:,nthComp),LDA_15score_Store,'UniformOutput',false);
score15=[score15;LDA_Comp];
end
score15=score15(finalSelect,:);
score15(pFinal==0)={[]};
score15=fillEmptyAsNaN(score15);
designMatrixCP = cell(5,15,3); % 5 ani 15 d 3 CP


WAPU_valid =(whichAnimalPerUnit(finalSelect));


for ani = 1:5
    ani
     temp=score15(sum(pFinal,2)>=greaterThanNd &  WAPU_valid==ani & umapWorstDistanceCri,:);
    for d = 1:15
        for CP = 1:max(nthComp)
        designMatrixCP{ani,d,CP}= [cell2mat(cellfun(@(x) mean(x(:,1:sampleEveryNTrial:end,CP),2)',temp(:,d),'Uni',0))';cell2mat(cellfun(@(x) mean(x(:,2:sampleEveryNTrial:end,CP),2)',temp(:,d),'Uni',0))'];
        end
    end
end


temp=designMatrixCP(:,:,1); % 
for ani=1:5
    for d=1:15
%
    temp{ani,d}=[designMatrixCP{ani,d,1} ;designMatrixCP{ani,d,2} ;designMatrixCP{ani,d,3}];
    temp{ani,d}= reshape(temp{ani,d},size(designMatrixCP{1},1),[]);
    end
end
designMatrixCP=temp;