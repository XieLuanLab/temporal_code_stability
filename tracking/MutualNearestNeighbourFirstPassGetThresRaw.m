function [FeatureThres,FeatDist,DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,XCorrPerTopChDistMatrix,FirLeadSecBy10XSample] = MutualNearestNeighbourFirstPassGetThresRaw...
    (DataForFeature,P2P,dayNumForEachUnit,MaskTrack_FarLoc,XCorrAvgPerTopChDistThres,FeatAndDist,keyPercentile)
% generate a 3 dimension distanceMatrix and return it,potentially inheritable by later programs. 
% auto discover thresholds on features looks at worst n%  threshold on input feat and input dist type for each distance informed MNN,
%% XCorrAvgPerChDistMatrix (mask thres)
numUnit = numel(DataForFeature);
numCh = size(DataForFeature{1},1);
NtopCh = min(4,numCh);

peakChList = zeros(numUnit,NtopCh);
for i = 1:numUnit
[~,I]=sort(P2P{i},'desc');
peakChList(i,:)=I(1:NtopCh);
end

X_flat = cellfun(@(x) reshape(x',numel(x(:)),1)',DataForFeature,'UniformOutput',0);
X_flat = double(cell2mat(X_flat));

XCorrDistFeatForTreeBuild      = cell(2,6); %dim1 [ raw or 1st diff] dim2 [ 5 common distance type(after top Ch alignment), and no alignment Xcorr]
XCorrPerTopChDistMatrix = zeros(numUnit);% this is used for mask thresholding. 

M11=nan(numUnit);M21=nan(numUnit);
M12=nan(numUnit);M22=nan(numUnit);
M13=nan(numUnit);M23=nan(numUnit);
M14=nan(numUnit);M24=nan(numUnit);
M15=nan(numUnit);M25=nan(numUnit);
M16=nan(numUnit);M26=nan(numUnit);

FirLeadSecBy10XSample = nan(numUnit);
parfor i = 1:numUnit
    i
    for j = 1:numUnit
        if j>i
        if MaskTrack_FarLoc(i,j)
            pooledChList = union(peakChList(i,P2P{j}(peakChList(i,:))~=0),peakChList(j,P2P{i}(peakChList(j,:))~=0));
            corrAllCh = zeros(1,numel(pooledChList));
            data3=DataForFeature{i};
            data4=DataForFeature{j};
            data1=data3(pooledChList,:);
            data2=data4(pooledChList,:);
            for ch = 1:numel(corrAllCh)
                x=data1(ch,:);
                y=data2(ch,:);
                corrAllCh(1,ch) = 1-max(xcorr(zscore(x'),zscore(y'),'normalized'));
            end
            meanCorrAllCh = mean(corrAllCh);
%% XCorrAvgPerTopChDistMatrix for mask Thres
                XCorrPerTopChDistMatrix(i,j)=meanCorrAllCh;            
  if meanCorrAllCh<=XCorrAvgPerTopChDistThres

%% XcorrDistXFlat; obtain 2 trees from xflat xcorr , xflatDiff xcorr
X_flat1=X_flat(i,:);
X_flat2=X_flat(j,:);
X_temp=X_flat2;
X_flat2(X_flat1==0 | X_flat2==0)=[];
X_flat1(X_flat1==0 | X_temp==0)=[];

x = zscore(X_flat1)';
y = zscore(X_flat2)';
M16(i,j)=1-max(xcorr(x,y,'normalized')); % counter ch loss with non-zero ele
x = zscore(diff(X_flat1))';
y = zscore(diff(X_flat2))';
M26(i,j)=1-max(xcorr(x,y,'normalized')); % counter ch loss with non-zero ele
%% Top N Ch Sample Alignment (with help from XCORR) 
upSampledShiftMax = 25
data1=resample(double(data1'),10,1)';
data2=resample(double(data2'),10,1)';
data3=resample(double(data3'),10,1)';
nT = size(data1,2);
r=xcorr2(data1,data2);
centerR = r(size(data1,1),:);
firstWaveLeadSecondWaveBy = find(centerR==max(centerR))-size(data1,2);
firstWaveLeadSecondWaveBy = firstWaveLeadSecondWaveBy(1);
FirLeadSecBy10XSample(i,j)=firstWaveLeadSecondWaveBy; 
if firstWaveLeadSecondWaveBy>upSampledShiftMax
    firstWaveLeadSecondWaveBy=upSampledShiftMax;
elseif firstWaveLeadSecondWaveBy<-upSampledShiftMax
    firstWaveLeadSecondWaveBy=-upSampledShiftMax;
end
data3=interp1(data3',1+firstWaveLeadSecondWaveBy:nT+firstWaveLeadSecondWaveBy,'linear','extrap')';
data3=interp1(0:nT-1,data3',0:10:nT-10,'linear','extrap')';

data3=data3';
data4=data4';
data3=single(data3(:)');
data4=data4(:)';
X_temp=data4;
data4(data4==0 | data3==0)=[];
data3(data3==0 | X_temp==0)=[];

                M11(i,j)=pdist2(data3,data4,'euclidean');
                M12(i,j)=pdist2(data3,data4,'cosine');
                M13(i,j)=pdist2(data3,data4,'cityblock');
                M14(i,j)=pdist2(data3,data4,'chebychev');
                M15(i,j)=pdist2(data3,data4,'correlation');    
data3=diff(data3);
data4=diff(data4);
                M21(i,j)=pdist2(data3,data4,'euclidean');
                M22(i,j)=pdist2(data3,data4,'cosine');
                M23(i,j)=pdist2(data3,data4,'cityblock');
                M24(i,j)=pdist2(data3,data4,'chebychev');
                M25(i,j)=pdist2(data3,data4,'correlation');    
                
            end
        end
        end
    end
end

XCorrDistFeatForTreeBuild{1,1}=M11; 
XCorrDistFeatForTreeBuild{1,2}=M12; 
XCorrDistFeatForTreeBuild{1,3}=M13; 
XCorrDistFeatForTreeBuild{1,4}=M14; 
XCorrDistFeatForTreeBuild{1,5}=M15; 
XCorrDistFeatForTreeBuild{1,6}=M16; 
XCorrDistFeatForTreeBuild{2,1}=M21; 
XCorrDistFeatForTreeBuild{2,2}=M22; 
XCorrDistFeatForTreeBuild{2,3}=M23; 
XCorrDistFeatForTreeBuild{2,4}=M24; 
XCorrDistFeatForTreeBuild{2,5}=M25; 
XCorrDistFeatForTreeBuild{2,6}=M26;
M11=[];M12=[];M13=[];M14=[];M15=[];M16=[];
M21=[];M22=[];M23=[];M24=[];M25=[];M26=[];
for i=1:2
    for j = 1:6
    MM = XCorrDistFeatForTreeBuild{i,j};
    MMisnanInitial = isnan(MM);
    MMisnanFinal = isnan(MM)&isnan(MM');
    MM(MMisnanInitial)=0;
    MM=MM+MM';
    MM(MMisnanFinal)=max(MM(:));

    XCorrDistFeatForTreeBuild{i,j}=MM.*(ones(size(MM))-eye(size(MM)));
    end
end
    MM = FirLeadSecBy10XSample;
    MMisnanInitial = isnan(MM);
    MMisnanFinal = isnan(MM)&isnan(MM');
    MM(MMisnanInitial)=0;
    MM=MM-MM';
    MM(MMisnanFinal)=nan;
    
FirLeadSecBy10XSample=MM;
XCorrPerTopChDistMatrix =XCorrPerTopChDistMatrix +XCorrPerTopChDistMatrix';
%% Standard Measures
MaskTrack_FarLocPerTopChXCorr= MaskTrack_FarLoc.*(XCorrPerTopChDistMatrix<=XCorrAvgPerTopChDistThres);
% input =  flattened waveform  Xnew
% output = FeatureThres
DistFeatForTreeBuild = cell(2,5,3); % so you have 3 dimensional distances to describe each tree 
FeatDist = cell(numel(DistFeatForTreeBuild)+numel(XCorrDistFeatForTreeBuild),size(FeatAndDist,1));
% dimension 1  = [ raw waveform , first derivative ] 
% dimension 2 =  [distance type used ]
% dimension 3 =  [type of dimension reduction used(None),PCA,UMAP]
numUnit = numel(dayNumForEachUnit);
mnn=1;
for wav = 1:2 % dimension 1 
 if wav==1 
  wave = X_flat;
 else
  wave = diff(X_flat,1,2);
 end
 for interest = 1:5 % dimension 2 
 switch interest   
          case 1
              strd = 'euclidean';%@distfunEUCL;% failed imput
          case 2 
              strd = 'cosine';%@distfunCos;
          case 3
              strd = 'cityblock';%@distfunCity;
          case 4
              strd = 'chebychev';%@distfunChev;
          case 5
              strd = 'correlation';%@distfunCORR;             
          otherwise
      end
 distMatrix = zeros(numUnit,numUnit,3);
 distMatrix(:,:,1) = pdist2(wave,wave,strd); % partial correction? 

 [~,score,~] = pca(wave(:,sum(wave==0,1)==0)); % have to use common chs 
 distMatrix(:,:,2) = pdist2(score(:,1:min(12,size(score,2))),score(:,1:min(12,size(score,2))),strd);

 [score, ~, ~, ~]=run_umap(wave(:,sum(wave==0,1)==0),'n_neighbors',numel(unique(dayNumForEachUnit))*6,'n_components',8,'min_dist',0.025,'metric',strd,'cluster_output','none','verbose','none');
 
 distMatrix(:,:,3) = pdist2(score,score);

for combine_dist = 1:3 % dimension 3
distCurrent  = distMatrix(:,:,combine_dist);
distCurrent(MaskTrack_FarLocPerTopChXCorr==0)=max(distCurrent(:));
DistFeatForTreeBuild{wav,interest,combine_dist}=distCurrent;
cluResult=MutualNearestNeighbour(distCurrent,true(size(distCurrent)),dayNumForEachUnit);

for feat = 1:size(FeatAndDist,1)
    FeatValues = FeatAndDist{feat,1};
    if isempty(FeatValues)
    FeatValues = X_flat;
    end
    DistType   = FeatAndDist{feat,2};
    FeatDist{mnn,feat} = StepMaxAvgDistMeanStdPrc9095(FeatValues,DistType,cluResult,keyPercentile);
end
    mnn = mnn +1;

end
end
end


for i=1:2
    for j =1:6

distCurrent  = XCorrDistFeatForTreeBuild{i,j}; % this dist not regulated by MaskTrack_FarLocPerTopChXCorr because it has already by implicitly incoporatexdin calculation steps above.
cluResult=MutualNearestNeighbour(distCurrent,true(size(distCurrent)),dayNumForEachUnit);

for feat = 1:size(FeatAndDist,1)
    FeatValues = FeatAndDist{feat,1};
    if isempty(FeatValues)
    FeatValues = X_flat;
    end
    DistType   = FeatAndDist{feat,2};
    FeatDist{mnn,feat} = StepMaxAvgDistMeanStdPrc9095(FeatValues,DistType,cluResult,keyPercentile);
end
    mnn = mnn +1;
    
    end
end



FeatureThres = cellfun(@(x) x(3,1),FeatDist); % the extracted 3 correspond to the dist4by3(3,:) = prctile(distNby3,keyPercentile,1); and 1 means mean within nearby day of same cluster and 3 means key percentile across clusters.
FeatureThres = mean(FeatureThres,1); % avg across all distances enumerated above 
end









function dist4by3 = StepMaxAvgDistMeanStdPrc9095(FeatValues,DistType,clu_Label,keyPercentile)
[NUPC,WhichClu]=histcounts(clu_Label,'BinMethod','integers');
validClu = WhichClu(NUPC>=2)+0.5;
distNby3 = zeros(numel(validClu),3);

distMatrix = pdist2(FeatValues,FeatValues,DistType);
for i = 1:numel(validClu)
subMatrix = distMatrix(clu_Label==validClu(i),clu_Label==validClu(i));

distNby3(i,1) = mean(diag(subMatrix,1));
distNby3(i,2) = mean(subMatrix(:));
distNby3(i,3) = max(subMatrix(:)); 

end
dist4by3 = zeros(4,3);
dist4by3(1,:) = mean(distNby3,1);
dist4by3(2,:) = std(distNby3,[],1);
dist4by3(3,:) = prctile(distNby3,keyPercentile,1); % 92.5% is used for across days
dist4by3(4,:) = prctile(distNby3,95,1);
end

