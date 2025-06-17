function [rawZ,maskMatrixHard,FeatureThres,maskMatrixLoose]      = MutualNearestNeighbour2ndWithThresRaw...
    (DataForFeature,DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,dayNumForEachUnit,MaskTrack_FarLocPerTopChXCorr,FeatAndDist,thres)
X_flat = cellfun(@(x) reshape(x',numel(x(:)),1)',DataForFeature,'UniformOutput',0);
X_flat = double(cell2mat(X_flat));
%% calculate the Joint FeatThres As maskMatrixHard to be taken in to account in the MNN funcitno 
numUnit = numel(dayNumForEachUnit);
distMatrix = zeros(numUnit,numUnit,size(FeatAndDist,1));
for feat = 1:size(FeatAndDist,1)
    FeatValues = FeatAndDist{feat,1};
    if isempty(FeatValues)
    FeatValues = X_flat;
    end
    DistType   = FeatAndDist{feat,2};
    distMatrix(:,:,feat) = pdist2(FeatValues,FeatValues,DistType);
    distMatrix(:,:,feat) = distMatrix(:,:,feat)<thres(feat);
end
maskMatrixHard = prod(distMatrix,3) & MaskTrack_FarLocPerTopChXCorr;
rawZ =[];
%% generate a 4 dimension distanceMatrix form MNN intialized trees , get 2sigma thresholds 
% input =  flattened waveform  Xnew
% output = FeatureThres
FeatDist = cell(numel(DistFeatForTreeBuild)+numel(XCorrDistFeatForTreeBuild),size(FeatAndDist,1));
% dimension 1  = [ raw waveform , first derivative ] 
% dimension 2 =  [distance type used ]
% dimension 3 =  [type of feature used , PCA]
% dimension 4(5) = [type of link used ] % since dimension 4 is a singleton dimension. not


mnn=1;
for wav = 1:2 % dimension 1 
for interest = 1:5 % dimension 2 
for combine_dist = 1:3 % dimension 3
distCurrent  = DistFeatForTreeBuild{wav,interest,combine_dist};
clu_Label=MutualNearestNeighbour(distCurrent,maskMatrixHard,dayNumForEachUnit);
[Z,~] = Link4Dist3onDistanceMatrix(distCurrent,clu_Label);
rawZ = [rawZ;Z];
for feat = 1:size(FeatAndDist,1)
    FeatValues = FeatAndDist{feat,1};
    if isempty(FeatValues)
    FeatValues = X_flat;
    end
    DistType   = FeatAndDist{feat,2};
    FeatDist{mnn,feat} = StepMaxAvgDistMeanStdPrc9095(FeatValues,DistType,clu_Label);
end
mnn = mnn +1;
end
end
end

for i=1:2
    for j =1:6
distCurrent  = XCorrDistFeatForTreeBuild{i,j}; % this dist not regulated by MaskTrack_FarLocPerTopChXCorr because it has already by implicitly incoporatexdin calculation steps above.
clu_Label=MutualNearestNeighbour(distCurrent,maskMatrixHard,dayNumForEachUnit);
[Z,~] = Link4Dist3onDistanceMatrix(distCurrent,clu_Label);
rawZ = [rawZ;Z];
for feat = 1:size(FeatAndDist,1)
    FeatValues = FeatAndDist{feat,1};
    if isempty(FeatValues)
    FeatValues = X_flat;
    end
    DistType   = FeatAndDist{feat,2};
    FeatDist{mnn,feat} = StepMaxAvgDistMeanStdPrc9095(FeatValues,DistType,clu_Label);
end
mnn = mnn +1;
    
    end
end





FeatureThres = cellfun(@(x) x(4,3),FeatDist);% 4,3 means max across cluster max within cluster. max within clu (across day) estimate expanded error. 
% Across All MNN 
FeatureThres = max(FeatureThres,[],1);      % a further max here mean max across distances (42) 
FeatureThres(FeatureThres<thres)=thres((FeatureThres<thres));
distMatrix = zeros(numUnit,numUnit,size(FeatAndDist,1));
for feat = 1:size(FeatAndDist,1)
    FeatValues = FeatAndDist{feat,1};
    if isempty(FeatValues)
    FeatValues = X_flat;
    end
    DistType   = FeatAndDist{feat,2};
    distMatrix(:,:,feat) = pdist2(FeatValues,FeatValues,DistType);
    distMatrix(:,:,feat) = distMatrix(:,:,feat)<FeatureThres(feat);
end
maskMatrixLoose = prod(distMatrix,3);
end

function dist4by3 = StepMaxAvgDistMeanStdPrc9095(FeatValues,DistType,clu_Label)
[NUPC,WhichClu]=histcounts(clu_Label,'BinMethod','integers');
validClu = WhichClu(NUPC>=2)+0.5;
distNby3 = zeros(numel(validClu),3);

distMatrix = pdist2(FeatValues,FeatValues,DistType);
for i = 1:numel(validClu)
subMatrix = distMatrix(clu_Label==validClu(i),clu_Label==validClu(i));
% Within Cluster 
numDiag=numel(diag(subMatrix,1));
distNby3(i,1) = sum(diag(subMatrix,1))/numDiag;
distNby3(i,2) = mean(subMatrix(:));
distNby3(i,3) = max(subMatrix(:)); 

end
dist4by3= zeros(4,3);
% Across Cluster
if numel(validClu)>=1
dist4by3(1,:) = mean(distNby3,1);
dist4by3(2,:) = std(distNby3,[],1);
dist4by3(3,:) = prctile(distNby3,90,1);
dist4by3(4,:) = max(distNby3,[],1);
end
end
