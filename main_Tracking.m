
addpath(genpath(pwd))
javaaddpath([pwd filesep 'umapAndEpp' filesep 'util' filesep 'suh.jar']);
%% Main
cd('shank2')
clear all; close all; % 
skipOneFile = 0; %used for DenoisePhoto()
trueFs = 25; % used for DenoisePhoto()
DenoisePhoto() % 
mkdir('results')
save('results/remain.mat','remain_index');
global ccg_PerSession withinSesValidCellId
ccg_PerSession = CCG_PerSession; % to be used when a button named CCG is pressed in the main plot.
withinSesValidCellId=WithinSesValidCellId;
WithinSesValidCellId=[]
CCG_PerSession=[];
P2PForFeature=P2PForFeature(remain_index==1);
DataForFeature = DataForFeature(remain_index==1);
St_WaveSeg = cellfun(@(x) x(:,10:59),DataForFeature,'UniformOutput',0); % 10:59  for Fs =25 K, spike peak at 31th pt.

maxP2P = double(cellfun(@max,P2PForFeature));

dayNumForEachUnit = dayNumForEachUnit(remain_index==1);
FWHM_summary      = FWHM_summary(remain_index==1);
crossP            = double(crossP(remain_index==1,:));
acg_metrics.acg_narrow=acg_metrics.acg_narrow(:,remain_index==1);
acg_metrics.acg_wide=acg_metrics.acg_wide(:,remain_index==1);
PeakChList= PeakChList(remain_index==1);
PosNegUnit = PosNegUnit(remain_index==1);
FRForFeature = FRForFeature(remain_index==1);

FeatAndDist = cell(4,2);
FeatAndDist{1,1} =maxP2P;       FeatAndDist{1,2} = 'euclidean';%
FeatAndDist{2,1} =FWHM_summary; FeatAndDist{2,2} = 'euclidean';
FeatAndDist{3,1} =[];           FeatAndDist{3,2} = 'corr';%
FeatAndDist{4,1} =crossP;       FeatAndDist{4,2} = @distfunXPPS;  % 
FeatAndDist{5,1} =crossP;       FeatAndDist{5,2} = @distfunX3PPS; %  
%% this assumes 2x16 design, and a wrongly resized channel map (historical reasons) must pause here and work on other electrode design
Ch_Map  = [result.Ch_Map_2(1:2,:) result.Ch_Map_2(3:4,:)];

debt_Ch_Map = Ch_Map;
debt_Ch_Map(debt_Ch_Map==0)=1;
options.debt_Ch_Map=debt_Ch_Map;
options.Ch_Map = Ch_Map;
WC10 = cellfun(@(x) regionprops(true(size(options.Ch_Map)), (x(options.debt_Ch_Map).^10).*(options.Ch_Map>0), 'WeightedCentroid'),P2PForFeature);
WC10 = cat(1,WC10.WeightedCentroid);
WC2= cellfun(@(x) regionprops(true(size(options.Ch_Map)), (x(options.debt_Ch_Map).^2).*(options.Ch_Map>0), 'WeightedCentroid'),P2PForFeature);
WC2 = cat(1,WC2.WeightedCentroid);


distThresEleStrongest= 0.75;        % only apply at 2nd itr
distThresEleAdjThre= 1.5; %         apply to nearby ses
distThresEleStrong= 2.25; %on average amount of drift from best knowledge.
distThresEleWeak = 4.5; %absolutely unnecessary to proceed if the drift is that much.  when respected, this save processing time along the way.
XCorrAvgPerTopChThresStrongest = 0.1; % only apply at 2nd itr
XCorrAvgPerTopChThresAdj = 0.2;    % apply to those at nearby ses
XCorrAvgPerTopChThresWeak = 0.3;      % absolutely respected thres
XPPSThresStrongest = 0.6;             % only apply at 2nd itr


MaskTrack_FarLoc =   (pdist2(WC2,WC2)<distThresEleWeak) & (pdist2(WC10,WC10)<distThresEleWeak)...
    &                ((pdist2(WC2,WC2)<distThresEleStrong) | (pdist2(WC10,WC10)<distThresEleStrong));
MaskTrack_NearLoc = (pdist2(WC2,WC2)<distThresEleWeak) & (pdist2(WC10,WC10)<distThresEleWeak)...
    &                ((pdist2(WC2,WC2)<distThresEleStrongest) | (pdist2(WC10,WC10)<distThresEleStrongest));
MaskTrack_AdjSesLoc = (pdist2(WC2,WC2)>=distThresEleAdjThre) & (pdist2(WC10,WC10)>=distThresEleAdjThre) & (pdist2(dayNumForEachUnit,dayNumForEachUnit)<=1);
tic
[FeatureThres,FeatDist,DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,XCorrPerTopChDistMatrix,FirLeadSecBy10XSample ] = MutualNearestNeighbourFirstPassGetThresRaw...
    (St_WaveSeg,P2PForFeature,dayNumForEachUnit,MaskTrack_FarLoc,XCorrAvgPerTopChThresWeak,FeatAndDist,95); % auto discover thresholds on features through building many trees and looks at worst n%
toc
% delete(gcp)
MaskTrack_AdjSesXCorr = (XCorrPerTopChDistMatrix>=XCorrAvgPerTopChThresAdj)  &  (pdist2(dayNumForEachUnit,dayNumForEachUnit)<=1) ;
MaskTrack_FarLoc(MaskTrack_AdjSesLoc)=0;MaskTrack_FarLoc(MaskTrack_AdjSesXCorr)=0;
MaskTrack_NearLoc(MaskTrack_AdjSesLoc)=0;MaskTrack_NearLoc(MaskTrack_AdjSesXCorr)=0;
XPPS  = pdist2(crossP,crossP,@distfunXPPS);
MaskTrack_XPPS    = XPPS < XPPSThresStrongest;
MaskTrack_FarLocPerTopChXCorrWeakWadj= MaskTrack_FarLoc.*(XCorrPerTopChDistMatrix<XCorrAvgPerTopChThresWeak);
MaskTrack_NearLocPerTopChXCorrStrongXPSSWadj= MaskTrack_NearLoc.*(XCorrPerTopChDistMatrix<XCorrAvgPerTopChThresStrongest).*MaskTrack_XPPS;

% important, this FeatureThres is meaningless or harmful when
% when total number of units/matched units are small
% which means if certain shanks have close to none recording, you should
% ignore it completely. The first and second feat also vulnerable to change
% in available chs, assign reasonable value (overwrite) your self.
if sum(isnan(FeatureThres))>0
    FeatureThres=[10 1.5 0.4 1 1];
end
[rawZ,maskMatrixHard,FeatureThresLoose,maskMatrixLoose]      = MutualNearestNeighbour2ndWithThresRaw(St_WaveSeg,DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,dayNumForEachUnit,MaskTrack_FarLocPerTopChXCorrWeakWadj,FeatAndDist,FeatureThres);
% new MNN with past thres, and estimate a loose thres that consider error
% exansion across time, form MNN initalized trees.
compositionForest = CompositonForestCreation(rawZ); %% a third pass on MNN result but without MNN
% see how high can you go up each tree with the loose thres
FinalPurityForest   = getFinalPurity(maskMatrixLoose,compositionForest);
AllPartition        = getPartition(FinalPurityForest,compositionForest);
% now how about re-initialize the tree , use autothreshold generatered
% clusters to initialize
[rawZthres] = ThresholdForestRaw(DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,AllPartition);

options.colorList=[0 0 1];
options.Yoff= 200;
options.Xoff= 60;
options.Visible = 'on';
options.subplotRow = 3;
options.subplotCol = 2;
options.subplotContent = cell(options.subplotRow,options.subplotCol); % specify subplot content
options.dayNumForEachUnit=dayNumForEachUnit;
options.requireParentNodeContinuousInTime=0; 

maxXPPSwhenPeakChDiffer=0;
preNumClu  = numel(St_WaveSeg);
numClu=preNumClu-1;
itr = 1;
outputArray = cell(1);
outputArray{1} = logical(eye(numel(St_WaveSeg)));
search =  2 ; % search 1 means a default situation where everything is unknown where search 2 incoporate priors.
%% usually worth a uisave at this point
save('matlab.mat')
%% Initial Start (skip next section)
curPathTemp=pwd;
shankStr = [curPathTemp(end-5:end) '/']  ; mkdir('videoFol')
cd ..
options.usePromptHistory = 1; % usage unclear
global PromptHistory; PromptHistory.prompt=cell(1,1);PromptHistory.composition=cell(1,1);PromptHistory.curPos =2;
PromptHistory.StartingCurPosItr=zeros(1,10);

%% Start from saved files maybe (only run one of the following: this section or the previous seciton )
% load('matlab.mat')
% options.dayNumForEachUnit=dayNumForEachUnit;% due to some past mistakes.
% curPathTemp=pwd;
% shankStr = [curPathTemp(end-5:end) '/']  ; mkdir('videoFol')
% cd ..
% addpath(genpath(pwd))
% options.usePromptHistory = 1;
% global PromptHistory; PromptHistory.prompt=cell(1,1);PromptHistory.composition=cell(1,1);PromptHistory.curPos =2;
% PromptHistory.StartingCurPosItr=zeros(1,10);

%% Iteration Loop
while numClu<preNumClu

    preNumClu=numClu;
    preUnitCount = size(outputArray{itr},1);

    outputArrayRecrur = outputArray{itr};

    WaveSeg=St_WaveSeg;
    P2P = P2PForFeature;


    WC = cellfun(@(x) regionprops(true(size(Ch_Map)), (x(debt_Ch_Map).^10).*(Ch_Map>0), 'WeightedCentroid'),P2P);
    WC = cat(1,WC.WeightedCentroid);
    PromptHistory.StartingCurPosItr(itr)=PromptHistory.curPos;

    itr = itr +1;


    numUnit = numel(WaveSeg);
    DataForFeatureAvgPer32days = WaveSeg;
    if itr==2
        plot_AcrossDay_units(DataForFeatureAvgPer32days,acg_metrics,dayNumForEachUnit,FRForFeature,options,shankStr)
    end
    currentLabel = 1:numel(DataForFeatureAvgPer32days);
    outputClusterSet=eye(numel(DataForFeatureAvgPer32days));
    cycle_recursive_result=0; % temporarily stores number of clusters returned for each forest
    ami_recursive_result=0;  %  temporarily stores AMI value for each forest
    runTime_result = 0;      % temporarily stores run time value for each forest
    promptResult = cell(1,1); % temporarily stores number of prompt user has to provide to the algorithm for each forest

    % Flatten waveforms , generate original forest

    corrMatrix = DistFeatForTreeBuild{1,5,1};
    if itr>2 % extract user input for easier program understanding.
        [MaskTrack_pastdecisions,badCompFromInput]= decodeUserPromptForMaskMatric(numUnit,PromptHistory,options);
    end
    % generate a forest of hierarchical clustering trees
    rng(1234)
    if itr==3
        MaskThisIter = MaskTrack_NearLocPerTopChXCorrStrongXPSSWadj;

        MaskThisIter(MaskTrack_pastdecisions==0)=0;
        [rawZ]= SystematicForestCreationRaw(DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,MaskThisIter);
        rawZstart = rawZ';
        options.defaultPrompt=-1;
        'Deafult 0 Expect to exercise Approval more'
        pause(1.5)
    elseif itr>3
        MaskThisIter = MaskTrack_FarLocPerTopChXCorrWeakWadj.*maskNN.*maskNoise;
        
        MaskThisIter(MaskTrack_pastdecisions==0)=0;
        XPPSinputThres=input(['max Within Clu XPPS when top Ch differ, same polarity is ' num2str(maxXPPSwhenPeakChDiffer) ' threshold at? : '])
        MaskThisIter((XPPS>=XPPSinputThres) & (pdist2(WC10,WC10)>0.5) & (pdist2(PosNegUnit,PosNegUnit)<0.5))=0;
        if exist('maskBadMerge')
            if input('ApplyMaskBadMerge:1/0')==1
                MaskThisIter(maskBadMerge==1)=-1;
            end
        end
        [rawZ]= SystematicForestCreationRaw(DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,MaskThisIter);
        rawZstart = rawZ';
        options.defaultPrompt=-1;
        'Deafult 0 Expect to exercise Approval more'
        pause(1.5)
    else
        MaskThisIter=maskMatrixLoose;
        MaskThisIter(MaskTrack_FarLocPerTopChXCorrWeakWadj==0)=0;
        rawZstart = rawZthres; % this varibale stores the 120 trees.
        options.defaultPrompt=-1;
        'Default 1 Expect to exercise Veto more'
        pause(1.5)
    end
    ranList = 1:max(size(rawZstart));% currently assuming all trees are used
    options.minNumFrameToAutoPlay=numel(unique(dayNumForEachUnit))*1.2;
    %% Perform random sampling
    expNum =1; % how many experiments to perform.
    treeN = numel(ranList); % how many trees to sample each experiment
    Final_numCluter = zeros(expNum,1); % assume 1 cycle required
    Final_AMI = zeros(expNum,1); % assume 1 cycle required
    Final_runTime=zeros(expNum,1);
    Final_PromptNum = zeros(expNum,1);
    Final_RandRecord = zeros(expNum,treeN);
    global figCells
    figCells=[];
   
    %% Important To Set queryGroundTruthMatrix as empty for program to run through manual label branch.
    queryGroundTruthMatrix=[]
    %% fix random seed for reproducibility
    figure('Position',[70         115        1851         854])
    drawnow
    for exp = 1:expNum % loop through all experiments

        tic
        rawZ=rawZstart(ranList,:);

        Final_RandRecord(exp,:)=ranList;

        if search == 1 && itr==2
            decisionStorage = -1*ones(numUnit,numUnit); % initialize as a N by N -1 array to store current decision.
            decisionStorage = decisionStorage + eye(size(decisionStorage))*2; % intialize 1 to diagonal.
        elseif search == 2 && itr==2
            decisionStorage = -1*ones(numUnit,numUnit); % initialize as a N by N -1 array to store current decision.
            decisionStorage = decisionStorage + eye(size(decisionStorage))*2; % intialize 1 to diagonal.
            decisionStorage(MaskThisIter==0)=0;
            'assignMatrixLosse'
        elseif search == 2 && itr==3
            decisionStorage = -1*ones(numUnit,numUnit); % initialize as a N by N -1 array to store current decision.
            decisionStorage = decisionStorage + eye(size(decisionStorage))*2; % intialize 1 to diagonal.
            decisionStorage(MaskThisIter==0)=0;
            decisionStorage(MaskTrack_pastdecisions==0)=0;
            decisionStorage(MaskTrack_pastdecisions==1)=1;
            'assignMaskTrack_NearLocPerTopChXCorrStrong'
        elseif search == 2 && itr>=4
            decisionStorage = -1*ones(numUnit,numUnit); % initialize as a N by N -1 array to store current decision.
            decisionStorage = decisionStorage + eye(size(decisionStorage))*2; % intialize 1 to diagonal.
            decisionStorage(MaskThisIter==0)=0;
            decisionStorage(MaskTrack_pastdecisions==0)=0;
            decisionStorage(MaskTrack_pastdecisions==1)=1;
            'assignMaskTrack_FarLocPerTopChXCorrWeak'
        end

        % convert original tree datastructure to that suitable for C-FAR
        compositionForest = CompositonForestCreation(rawZ);


        PF_depth = zeros(1,numel(compositionForest)); % find the depth of all trees inside forest, this function might be unused
        parfor tr = 1:numel(compositionForest)
            PF_depth(tr) = compositionForest{tr}.depth;
        end
        [sortedHeight,sortHeight]=sort(PF_depth);

        compositionForest=compositionForest(sortHeight);
        PF_depth=sortedHeight;
        % ----------------- convert from tree/cell structure to multi-dimensional--
        % array
        % originally in composition Forest it storesd unit identity # 1 # 5 # 84 something like that indicating cluster identity
        % in the new version compositionForest_bin, each item is a logical array of size = total
        % unit count, the zero or one at each bit indicated whether certain # unit
        % is in this cluster or not. In the latest compostionForest_binArray, it further convert Forest_bin into a multi-D array for faster computation

        compositionForest_bin = cellfun(@(x) x.treefun(@(Y) array2logical(Y,numUnit)), compositionForest,'UniformOutput',0); % create one hot encoding from original node compositions.
        compositionForest=[];
        purityForest = cellfun(@(x) tree(x,-1), compositionForest_bin,'UniformOutput',0);% reinitialize -1 as unknown

        for tr = 1:numel(purityForest)
            for i = purityForest{tr}.findleaves
                purityForest{tr}=purityForest{tr}.set(i,1);
                %initialize all leaf node to 1
            end
        end
        maxNodeNumber = max(cellfun(@(x) numel(x.Parent),compositionForest_bin));
        compositionForest_binArray = cellfun(@(x) PaddedCell2Mat(x.Node,maxNodeNumber),compositionForest_bin,'UniformOutput',0); % convert datatype from cell to multi-D array
        compositionForest_bin=[];
        compositionForest_binArray = reshape(cell2mat(compositionForest_binArray),size(compositionForest_binArray{1},1),size(compositionForest_binArray{1},2),numel(compositionForest_binArray));

        PFArray_parent = uint32(alignedCell2Mat(cellfun(@(x) x.Parent,purityForest,'UniformOutput',0),0)); % convert from original tree data structure to a vectorized representation for speed improvment
        PFArray_purity = int8(alignedCell2Mat(cellfun(@(x) cell2mat(x.Node),purityForest,'UniformOutput',0),0)); % same purpose as the previous line


        maxDepth = max(PF_depth);
        PFArray_ptr = zeros( size(PFArray_purity,1),maxDepth+1,numel(PF_depth),'uint16');
        parfor tr = 1:numel(PF_depth)
            tr
            temp = nan(size(PFArray_purity,1),maxDepth+1);
            for node = 1:numel(purityForest{tr}.Node)
                pathToRoot = purityForest{tr}.findpathUp(node, 1);
                temp(node,1:numel(pathToRoot))=pathToRoot;
            end

            PFArray_ptr(:,:,tr)=temp;
        end
        temp=[];

        Tree_by_UnitID_store_leafInd = zeros(numel(purityForest),numUnit); % precompute the index of all leaf nodes across all the trees in the forest for faster retrival of this info later
        leafForest  = cellfun(@(x) x.findleaves,purityForest,'UniformOutput',0);
        for ctr = 1:numel(purityForest)
            [r,c]=find(compositionForest_binArray(leafForest{ctr},:,ctr));
            [~,ci] = sort(c);
            Tree_by_UnitID_store_leafInd(ctr,:)=leafForest{ctr}(r(ci));
        end
        purityForest=[];

        totalUnitSet = array2logical(1:numUnit,numUnit);
        outputClusterSet = false(1,numUnit);

        node_by_TreeMatrix =(repmat(1:size(compositionForest_binArray,1),size(compositionForest_binArray,3),1))';
        if itr>2
            if ~isempty(badCompFromInput{1}) % incoporate past prompts, set their min-ext to zeros.
                for badC=1:numel(badCompFromInput)
                    badcomp = badCompFromInput{badC};
                    testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                    PFArray_purity(testM)=0;
                end
            end
        end

        purityForest_ori = PFArray_purity; % store a copy of purity tree for faster reinitialization



        cycle = 0;
        verdictCountStore = zeros(1,1);
        verdictCountType = zeros(1,1);
        %----- Algorithm 1 FindOneBloc (Alg1)  referenced in the paper
        % Input = {T1, T2, . . . , Tr}  <- compositionForest_binArray
        % Output =   one node among the trees which corresponds to one cluster best compatible with the pairwise feedbac <- outputClusterSet(end+1,:)

        while 1 % this loop check if there is still non-empty trees remaining in the forest(whether all units have been assigned a cluster label), if so,  break
            if (all(sum(compositionForest_binArray(:,:,1)')==0))% when all branches have been "chopped", break
                break
            end

            remainingUnit  = find(totalUnitSet - sum(outputClusterSet,1));% find out which units havn't been output yet.
            sprintf('remaining units %0.2f %%',numel(remainingUnit)/numUnit*100)
            % Alg1.1 Start with a definite pure set S (in this case, a leaf node)
            [~,leftMostUnit]=min(WC(remainingUnit,1));
            S = remainingUnit(leftMostUnit);

            pureNodeIndexPerTree  = Tree_by_UnitID_store_leafInd(:,S);

            while 1 % this loop iteratively find the largest pure node for the outter loop to cut

                unionOfS  =false(1,numUnit);
                % Alg1.7.2 S' = take union across unit composition across all pure nodes.
                for tr = 1:numel(PF_depth)
                    unionOfS = logical(unionOfS + compositionForest_binArray(pureNodeIndexPerTree(tr),:,tr));
                end

                % Alg 1.1 find minimal extensions for all the trees
                testM = squeeze(sum(compositionForest_binArray(:,unionOfS,:),2))==sum(unionOfS);
                testM = testM.*node_by_TreeMatrix;
                for tr = 1:numel(PF_depth)% loop through each tree to get minExtPerTree
                    testM(testM(:,tr)>=pureNodeIndexPerTree(tr),tr)=0;
                end

                MinExtNodeIndicesPerTree =max(testM);
                MinExtNodeIndicesPerTree(MinExtNodeIndicesPerTree==0)=1;
                % Alg 1.2 Determine the purity of min Extensions nodes.
                % make queries and update purityforest , update decisionstorage
                if itr==2
                    [minExtPurity,PFArray_purity,decisionStorage,verdictCount] = queryForestNodeChildren_binWithinDayFirstIterLtd(compositionForest_binArray,PFArray_purity,PFArray_parent,MinExtNodeIndicesPerTree,decisionStorage,options,queryGroundTruthMatrix,DataForFeatureAvgPer32days,shankStr,corrMatrix);
                else
                    [minExtPurity,PFArray_purity,decisionStorage,verdictCount] = queryForestNodeChildren_bin(compositionForest_binArray,PFArray_purity,PFArray_parent,MinExtNodeIndicesPerTree,decisionStorage,options,queryGroundTruthMatrix,DataForFeatureAvgPer32days,shankStr);
                end
                verdictCountStore(end+1)=verdictCount;
                verdictCountType(end+1) = 1;
                whichTree = find(minExtPurity==1);
                % Alg1.3 find out the value of m ( how many pure nodes )
                %             sum(minExtPurity)

                if any(minExtPurity==1 & MinExtNodeIndicesPerTree ==1) % if root node of any tree becomes pure
                    largestPureNode=1;
                    largestPureNodeComposition = find(compositionForest_binArray(largestPureNode,:,1));

                    decisionStorage(largestPureNodeComposition,largestPureNodeComposition)=1;
                    outputClusterSet(end+1,:)= compositionForest_binArray(largestPureNode,:,1);
                    break;
                end
                % Alg1.4 do a switch statement on m
                switch sum(minExtPurity)
                    case 0 % this is not described in Alg1, but implied, when there is no pure nodes among the extensions, then output the current S'
                        outputClusterSet(end+1,:) = unionOfS;
                        break;
                        % Alg1.5 do binary search along a single tree until the largest
                        % node B with purity of 1 is found
                    case 1
                        [largestPureNode,PFArray_purity,decisionStorage,verdictCount] =  queryForestNodeParent_bin(3,compositionForest_binArray,PFArray_purity,PFArray_parent,PFArray_ptr(:,:,whichTree),MinExtNodeIndicesPerTree,whichTree,decisionStorage,options,queryGroundTruthMatrix,DataForFeatureAvgPer32days,shankStr);
                        outputClusterSet(end+1,:)= compositionForest_binArray(largestPureNode,:,whichTree);
                        verdictCountStore(end+1) = verdictCount;
                        verdictCountType(end+1)  = 2;
                        break;
                    otherwise
                        % Alg1.8 Repeat with updated pure nodes
                        pureNodeIndexPerTree(whichTree) = MinExtNodeIndicesPerTree(whichTree); % update some pureNodeIndices with minExt node indices.
                end
            end
            cycle = cycle+1;
            ['loop exit: ' num2str(cycle)]
            compositionForest_binArray(:,outputClusterSet(end,:),:)=0;% pruneForest
            PFArray_purity=purityForest_ori; % re-intialize

        end
        %----------- Store result

        outputClusterSet(1,:)=[];
        outputArray{itr} = outputClusterSet;
        numNewUnit = size(outputClusterSet,1);

        if search==2 && itr ==2
            decisionStorage(maskMatrixLoose==0)=-1;
            decisionStorage(MaskTrack_FarLocPerTopChXCorrWeakWadj==0)=-1;
        elseif search == 2 && itr==3
            decisionStorage(MaskTrack_NearLocPerTopChXCorrStrongXPSSWadj==0)=-1;
        elseif search == 2 && itr>=4
            decisionStorage(MaskTrack_FarLocPerTopChXCorrWeakWadj==0)=-1;
        end


        promptResult{1} = verdictCountStore;

        cycle_recursive_result(1)=cycle;
             runTime_result(1) = 0;

        % accumulate results for this experiment
        Final_numCluter(exp,1:numel(cycle_recursive_result))=cycle_recursive_result;
       
        Final_runTime(exp,1:numel(runTime_result))=runTime_result;
        Final_PromptNum(exp,1:numel(cellfun(@sum,promptResult)))=cellfun(@sum,promptResult);

    end

    close all
    if itr>=3 % obtain maxXPPSwhenPeakChDiffer for user approved clusters
        for rep = 1:size(outputClusterSet,1)
            UnitsOfThatCluster = find(outputClusterSet(rep,:)==1);
            XPPSthisClu = XPPS(UnitsOfThatCluster ,UnitsOfThatCluster );
            WC10thisClu = pdist2(WC10(UnitsOfThatCluster,:),WC10(UnitsOfThatCluster,:));
            PosNegthisClu = pdist2(PosNegUnit(UnitsOfThatCluster,:),PosNegUnit(UnitsOfThatCluster,:));
            XPPSthisClu(WC10thisClu<0.5)=0; % we only care about a separate XPPS thres when peak ch differ.
            XPPSthisClu(PosNegthisClu~=0)=0; % we only care about a separate XPPS thres when two units have same polarity
            maxXPPSwhenPeakChDiffer=max(maxXPPSwhenPeakChDiffer,max(XPPSthisClu,[],'all'));
        end
    end
    if itr==3 % rule out some obvious "should not merge" cases so that later iterations are easier. first show singletons from past itr, ask if it's noise
        % for non-noise, do a 6NN (sum zscore of all ch waveform xcorr and WC2 distance)  search, no trees, just show 1 unit per cluster and their 6 nearest neighbours (respect later itr cutoffs)
        % users are expected to answer which of the 6 plots is acceptable (this decision however will not be respected, non acceptable ones will be respected), then we
        % update PromptHistory accordingly so that later iterations are easier.
        DistNN = zscore(pdist2(WC2,WC2),0,'all')+zscore(XCorrDistFeatForTreeBuild{1,6},0,'all');
        DistNN(MaskTrack_FarLocPerTopChXCorrWeakWadj==0)=NaN; % note this is next itr standard .
        [MaskTrack_pastdecisions,badCompFromInput]= decodeUserPromptForMaskMatric(numUnit,PromptHistory,options)
        DistNN(MaskTrack_pastdecisions==0)=NaN;
        RepUnitThatCluster = zeros(size(outputClusterSet,1),1)
        for rep = 1:numel(RepUnitThatCluster)
            UnitsOfThatCluster = find(outputClusterSet(rep,:)==1);
            RepUnitThatCluster(rep)=UnitsOfThatCluster(1);
        end
        DistNNrep = DistNN(RepUnitThatCluster,RepUnitThatCluster);
        for badC=1:numel(badCompFromInput)
            badC
            if ~isempty(badCompFromInput{badC})
                concernedClusters = find(sum(outputClusterSet(:,badCompFromInput{badC}),2))
                DistNNrep(concernedClusters,concernedClusters)=NaN;
            end
        end
        for rep = 1:numel(RepUnitThatCluster)
            DistNNrep(rep,rep)=NaN;
        end
        NumUnitThatClu = sum(outputClusterSet,2);
        PotNoiseList = find(NumUnitThatClu==1);
        NoiseList = []
        NoiseListIdx = []
        for potNoise = 1:numel(PotNoiseList)
            currentComposition=find(outputClusterSet(PotNoiseList(potNoise),:));
            [fig_handle] = videofig(shankStr,currentComposition,numel(currentComposition), @redrawAcrossDayUnits); % import
            redrawAcrossDayUnits(1,currentComposition,shankStr);
            drawnow;
            verdict = input([' Noise =? 1/yes, 0/no :']);
            if numel(verdict)==1 && sum(verdict==[0 1 ])==0
                verdict = input(['past input fault, was:' num2str(verdict) ' enter new one: ']);
            end
            if ~isempty(verdict) && verdict==1
                NoiseList(end+1)=currentComposition;
                NoiseListIdx(end+1)=potNoise;

                DistNNrep(PotNoiseList(potNoise),:)=NaN;
                DistNNrep(:,PotNoiseList(potNoise))=NaN;
            end
            try
                close(fig_handle)
                fig_handle=[];
            catch
            end
        end
        DistNNrep(PotNoiseList(NoiseListIdx),:)=NaN;
        DistNNrep(:,PotNoiseList(NoiseListIdx))=NaN;
        maskNoise = ones(size(DistNN));
        maskNoise(NoiseList,:)=0;
        maskNoise(:,NoiseList)=0;

        figure('Position',[70         115        1851         854])
        drawnow
        maskNN  = ones(size(DistNN));
        % DataForFeatureAvgPer32daysNorm=DataForFeatureAvgPer32days;
        % for DU=1:numel(DataForFeatureAvgPer32days)
        % DataForFeatureAvgPer32daysNorm{DU}=DataForFeatureAvgPer32days{DU}/crossP(DU,PeakChList(DU))*100;
        % end
        for rep = 1:numel(RepUnitThatCluster)
            rep
            while any(~isnan(DistNNrep(rep,:)))
                [V,sI]=sort(DistNNrep(rep,:),'desc','MissingPlacement','last');
                if sum(~isnan(DistNNrep(rep,:)))>6
                    finalsel = sI(1:6);
                else
                    finalsel = sI(~isnan(V));
                end
                options.currentComposition=RepUnitThatCluster(union(finalsel,rep));
                for subplotRow = 1:options.subplotRow
                    for subplotCol = 1:options.subplotCol
                        subplotId =(subplotRow-1)*options.subplotCol+subplotCol;
                        if subplotId>numel(finalsel)
                            continue
                        end
                        options.subplotContent{subplotRow,subplotCol}=[RepUnitThatCluster(rep) RepUnitThatCluster(finalsel(subplotId))];
                        subplot_tight(options.subplotRow,options.subplotCol,subplotId,[0.01 0.01]);
                        plot_two_units(subplotRow,subplotCol,options,DataForFeatureAvgPer32days)
                    end
                end
                drawnow
                verdict = input(['With Leap of Faith, Which plots might be same untis =? :']);
                badSuplots = setdiff(1:numel(finalsel),verdict);

                DistNNrepBk = DistNNrep;% you can regret for 1 step; Ctl C out .
                maskNNBk = maskNN;
                DistNNrep(rep,finalsel)=NaN;
                DistNNrep(finalsel,rep)=NaN;% this is path dependent.
                for badE = 1:numel(badSuplots)
                    AC=allcomb(find(outputClusterSet(rep,:)),find(outputClusterSet(finalsel(badSuplots(badE)),:)),'matlab');
                    maskNN((AC(:,1)-1)*numUnit+AC(:,2))=0;
                    maskNN((AC(:,2)-1)*numUnit+AC(:,1))=0;
                end
                clf
                if ~isempty(verdict) && ~isequal(verdict,0)
                    XPPSthisClu = XPPS(RepUnitThatCluster(rep) ,finalsel(verdict));
                    WC10thisClu = pdist2(WC10(RepUnitThatCluster(rep),:),WC10(finalsel(verdict),:));
                    PosNegthisClu = pdist2(PosNegUnit(RepUnitThatCluster(rep),:),PosNegUnit(finalsel(verdict),:));
                    XPPSthisClu(PosNegthisClu~=0)=0; % we care more about XPPS when Polarity is the same.
                    XPPSthisClu(WC10thisClu<0.5)=0; % we care more about XPPS when peak ch differ.
                    maxXPPSwhenPeakChDiffer=max(maxXPPSwhenPeakChDiffer,max(XPPSthisClu,[],'all'));
                    %     break
                end
            end
        end
        close all;

    end
    combinedResult = [Final_AMI(:,1) Final_numCluter(:,1) Final_PromptNum(:,1) toc];
    saveStr = char(datetime);
    save([shankStr 'results/' replace(saveStr,{' ','-',':'},'_')],'Final_RandRecord','numUnit','expNum'...
        ,'treeN','rawZstart','combinedResult','outputArray','decisionStorage','WaveSeg',...
        'dayNumForEachUnit','options','sortHeight','PromptHistory','MaskThisIter')
    if itr>=3
        save([shankStr 'results/' replace(saveStr,{' ','-',':'},'_')],'Final_RandRecord','numUnit','expNum'...
            ,'treeN','rawZstart','combinedResult','outputArray','decisionStorage','WaveSeg',...
            'dayNumForEachUnit','options','sortHeight','PromptHistory','MaskThisIter','NoiseList','maskNN','maxXPPSwhenPeakChDiffer')
    end
    numClu = Final_numCluter(:,1)
    'Ctrl+C at any time to stop going into next iteration'
end
%% Process results, load from relevant result file(first one)
cd([shankStr 'results'])
DataForFeatureAvgPer32days = WaveSeg;
ChMap = options.Ch_Map;
options=[];
Label = outputArray{end};
numClu = size(Label,1)
%% Plot overlaid result 0
maxSes = max(unique(dayNumForEachUnit));
options.colorList=winter(maxSes);
options.alpha= 0.6;
options.Yoff= 250;
options.Xoff= 50;
options.Visible = 'off';
options.subplotRow = 1;
options.subplotCol = 1;
options.subplotContent = cell(options.subplotRow,options.subplotCol); % specify subplot content
options.dayNumForEachUnit=dayNumForEachUnit;
options.Ch_Map = ChMap;
try
    rmdir('testplot2m','s')
catch
end
mkdir('testplot2m')
parfor i = 1:size(Label,1)
    i
    uList = Label(i,:);
    if sum(uList)>0
        figure('Position',[155          88        1108         599],'Visible',options.Visible)
        plot_two_units_simple(DataForFeatureAvgPer32days,find(uList),options)
        colormap(options.colorList)
        if sum(uList(NoiseList))>0
            colormap autumn
        end
        c=colorbar;
        c.Ticks = linspace(0, 1, maxSes) ; %Create 8 ticks from zero to 1
        c.TickLabels = num2cell(1:maxSes) ;    %R
        c.FontSize=14;
        title(sprintf('%d waveforms appeared on %d sessions',sum(uList),numel(unique(dayNumForEachUnit(uList)))))

        saveas(gcf,['testplot2m' filesep num2str(i) '.png'])
        close(gcf)
    end
end

