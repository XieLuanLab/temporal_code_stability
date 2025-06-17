function [minExtPurity,PFArray_purityOri,decisionStorage,verdictCount,BadCompInputStore] = queryForestNodeChildren_binWithinDayFirstIterLtd(compositionForest_binArrayOri,PFArray_purityOri,PFArray_parentOri,MinExtNodeIndicesPerTree,decisionStorage,options,GTmode,DataForFeatureAvgPer32days,shankStr,corrMatrix)
%                                                                                                                                           compositionForest_binArray,PFArray_purity,PFArray_parent,MinExtNodeIndicesPerTree,decisionStorage,options,queryGroundTruthMatrix
% This  function decide the purity of one node for each tree inside the
% whole forest.

% Input = composition forest, purity forest (parent and purity) , past
% decision, indices of node in question, possible ground truth matrix

% Output = purity of the nodes in question, updated purity forest, decision
% matrix, number of feedback required to get decision.

% This is a core process in Alg1 in Zhu, H., Li, X., Sun, L., He, F., Zhao, Z., Luan, L., ... & Xie, C. (2020). 
% Clustering with fast, automated and reproducible assessment applied to longitudinal neural tracking. arXiv preprint arXiv:2003.08533.

% This code also examplify the quote "Purity satisfies important inequalities
% imposed by set inclusion" by simultaneously
% eliminating "bad compositions" across all trees inside the forest. This
% is a key innovation for multi-tree setup. search "badcomp" for more
% information
% 'enter'
global PromptHistory

BadCompInputStore = [];
verdictCount = 0;
numTree = size(PFArray_purityOri,2);

clusterLabel=[]; % assume each minExt as a cluster, whose label is encoded in clusterLabel, 2D matrix, number of waveform-by-number of Tree
% # of cluster = # of trees, Then form super-cluster out of that and make a tree
for tr = 1:numTree % --1 loop through each tree
    clusterLabel(:,tr) = compositionForest_binArrayOri(MinExtNodeIndicesPerTree(tr),:,tr); % g
end
clusterLabel=clusterLabel==1;
testD = true(1,numTree);
for t = 1:numTree
    test = decisionStorage(clusterLabel(:,t),clusterLabel(:,t));
    testD(t) = ~any(test(:)==0);
end
if any(testD)
    clusterLabel=clusterLabel(:,testD);
    [UniqueclusterLabel,ia,ic] = unique(clusterLabel','rows');
    UniqueclusterLabel = UniqueclusterLabel';
    if size(UniqueclusterLabel,2)==1
        UniqueclusterLabel=repmat(UniqueclusterLabel,1,2);
    end
    
    
    %   size(UniqueclusterLabel,2)
    
    numUnit=size(corrMatrix,1);
    [Z, ~] = Link1Dist1onDistanceMatrixClusterLabelMatrix(corrMatrix,UniqueclusterLabel);
    compositionForest = CompositonForestCreation(Z);
    compositionForest{1}.Node=cellfun(@unique,compositionForest{1}.Node,'UniformOutput',false);
    compositionForest_bin = cellfun(@(x) x.treefun(@(Y) array2logical(Y,numUnit)), compositionForest,'UniformOutput',0); % create one hot encoding from original node compositions.
    compositionForest=[];
    purityForest = cellfun(@(x) tree(x,-1), compositionForest_bin,'UniformOutput',0);% reinitialize -1 as unknown
    %   for tr = 1:numel(purityForest)
    %       for i = purityForest{tr}.findleaves
    %           purityForest{tr}=purityForest{tr}.set(i,1);
    %           %initialize all leaf node to 1
    %       end
    %   end
    maxNodeNumber = max(cellfun(@(x) numel(x.Parent),compositionForest_bin));
    compositionForest_binArray = cellfun(@(x) PaddedCell2Mat(x.Node,maxNodeNumber),compositionForest_bin,'UniformOutput',0); % convert datatype from cell to multi-D array
    compositionForest_bin=[];
    compositionForest_binArray = reshape(cell2mat(compositionForest_binArray),size(compositionForest_binArray{1},1),size(compositionForest_binArray{1},2),numel(compositionForest_binArray));
    
    PFArray_parent = uint32(alignedCell2Mat(cellfun(@(x) x.Parent,purityForest,'UniformOutput',0),0)); % convert from original tree data structure to a vectorized representation for speed improvment
    PFArray_purity = int8(alignedCell2Mat(cellfun(@(x) cell2mat(x.Node),purityForest,'UniformOutput',0),0)); % same purpose as the previous line
    leaf=purityForest{1}.findleaves;
    purityForest{1}.Node(leaf)=[];
    purityForest{1}.Parent(leaf)=[];
    iterator = purityForest{1}.breadthfirstiterator;
    BFS_purity = -1*ones(1,numel(iterator));
    for node=1:numel(iterator)
        BFS_purity(node)=PFArray_purity(iterator(node));
        badCompFromInput=[];
        if BFS_purity(node)==-1
            currentComposition = find(compositionForest_binArray(iterator(node),:)); % get unit composition
            currentDecisionMatrix = decisionStorage(currentComposition,currentComposition); % obtain subMatrix from decision storage.
            numOfMinusOneEntry = sum(currentDecisionMatrix(:)==-1); % sum up all previously -1s
            try
                if options.requireParentNodeContinuousInTime==1
                    if any(diff(sort(options.dayNumForEachUnit(currentComposition)))>1)
                        currentDecisionMatrix(1)=0;
                    end
                end
            catch
            end
            if  numOfMinusOneEntry >=2 % decision Matrix of MinExtNode Compositon contains uncertainty, then intiate display query process
                %------------------------------------%obtain info and configure plotting options
                currentNodeChildren = find(PFArray_parent==iterator(node));% --2 identify the childern of the node in the current nodee
                leftChild = currentNodeChildren(1);
                rightChild = currentNodeChildren(2); % assume binary tree
                leftChildComp = sort(find(compositionForest_binArray(leftChild,:)));
                rightChildComp = sort(find(compositionForest_binArray(rightChild,:)));
                
                if numel(currentNodeChildren)>2 || numel(currentNodeChildren)==0
                    leftChild=iterator(node);
                    rightChild= iterator(node);
                    leftChildComp = currentComposition;
                    rightChildComp = currentComposition;
                end
                
                % might be better to assign -1
                options.currentComposition=currentComposition;
                options.subplotContent{1,1}=leftChildComp;
                options.subplotContent{1,2}=rightChildComp;
                options.subplotContent{2,1}=[leftChildComp(1) rightChildComp(1)];
                options.subplotContent{2,2}=[leftChildComp(end) rightChildComp(end)];
                options.subplotContent{3,1}=[leftChildComp(1) rightChildComp(end)];
                options.subplotContent{3,2}=[leftChildComp(end) rightChildComp(1)];
                
                if isempty(GTmode) % obtain feedback from users either by displaying two units 's waveform or by looking at ground truth
                    
                    if ~any(currentDecisionMatrix(:)==0)
                        
                        if options.defaultPrompt==1
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 14+numel(currentComposition)*0.04, 'TimerFcn', @pressOneEnter);
                        elseif options.defaultPrompt==2
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 1.5, 'TimerFcn', @pressTwoEnter);
                        elseif options.defaultPrompt==0
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 10, 'TimerFcn', @pressZeroEnter);
                        end
                        
                        
                        if options.usePromptHistory==1 && PromptHistory.curPos<=numel(PromptHistory.prompt)
                            
                            verdict = PromptHistory.prompt{PromptHistory.curPos};
                            PromptHistory.curPos=PromptHistory.curPos+1;
                        else
                            for subplotRow = 1:options.subplotRow
                                for subplotCol = 1:options.subplotCol
                                    subplot_tight(options.subplotRow,options.subplotCol,(subplotRow-1)*options.subplotCol+subplotCol,[0.01 0.01]);
                                    plot_two_units(subplotRow,subplotCol,options,DataForFeatureAvgPer32days)
                                    
                                end
                            end
                            displayChoiseButtons()
                            drawnow;
                            if sum(options.defaultPrompt==[0 1 2])~=0
                                start(t);
                                verdict = input(['combine=? 1/y, 0/n, 2/m, default ' num2str(options.defaultPrompt) ':']);
                                if (numel(verdict)==1 && sum(verdict==[0 1 2])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
                                    verdict = input(['past input fault, was wrong enter new one: ']);
                                end
                                PromptHistory.prompt{end+1}=verdict;
                                PromptHistory.composition{end+1}=currentComposition;
                                PromptHistory.curPos=PromptHistory.curPos+1;
                                currentComposition
                                verdict 
                                PromptHistory.curPos
                                stop(t);
                                delete(t);
                            else
                                verdict = input(['combine=? 1/y, 0/n, 2/m, default ' num2str(options.defaultPrompt) ':']);
                                if (numel(verdict)==1 && sum(verdict==[0 1 2])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
                                    verdict = input(['past input fault, was wrong enter new one: ']);
                                end
                                PromptHistory.prompt{end+1}=verdict;
                                PromptHistory.composition{end+1}=currentComposition;
                                PromptHistory.curPos=PromptHistory.curPos+1;
                                currentComposition
                                verdict 
                                PromptHistory.curPos
                            end
                        end
                        if isempty(verdict)
                            verdict = 0;
                        elseif numel(verdict)>1
                            verdict=real(verdict)+imag(verdict);
                            negativeElementsLocation = find(verdict(:,1)<0);
                            if ~isempty(negativeElementsLocation)
                                negativeElements = currentComposition(abs(verdict(negativeElementsLocation,1)));
                                restOfElements = setdiff(currentComposition,negativeElements);
                                badCompFromInput=[badCompFromInput;allcomb(negativeElements,restOfElements)];
                                verdict(negativeElementsLocation,:)=[];
                            end
                            if ~isempty(verdict)
                                badCompFromInput=[badCompFromInput;currentComposition(verdict)];
                            end
                            verdict  = 0;
                        end
                        
                        verdictCount = verdictCount+1;
                        
                        if verdict ==2
                            
                            if options.usePromptHistory==1 && PromptHistory.curPos<=numel(PromptHistory.prompt)
                                verdict = PromptHistory.prompt{PromptHistory.curPos};
                                PromptHistory.curPos=PromptHistory.curPos+1;
                            else
                                %                     plot_perDay_units(options);
                                uicontrol('Position',[800+950 750-600 60 48],'String','CCG',...
                                    'FontSize',14,'Callback', {@CCGgroupPlot,currentComposition,options.dayNumForEachUnit});%'BackgroundColor','r',
                                [fig_handle] = videofig(shankStr,currentComposition,numel(currentComposition), @redrawAcrossDayUnits); % import
                                redrawAcrossDayUnits(1,currentComposition,shankStr);
                                drawnow;
                                t = timer('ExecutionMode', 'singleShot', 'StartDelay', 0.4, 'TimerFcn', @pressEnter);
                                
                                if numel(currentComposition)<=options.minNumFrameToAutoPlay
                                    start(t);
                                    verdict = input(['combine=? 1/yes, 0/no :']);
                                    if (numel(verdict)==1 && sum(verdict==[0 1])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
                                        verdict = input(['past input fault, was wrong enter new one: ']);
                                    end
                                    stop(t);
                                else
                                    verdict = input(['combine=? 1/yes, 0/no :']);
                                    if (numel(verdict)==1 && sum(verdict==[0 1])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition)) 
                                        verdict = input(['past input fault, was wrong enter new one: ']);
                                    end
                                end
                                delete(t);
                                global figCells
                                if ~isempty(figCells) && ~isempty(figCells{1})
                                    for figC=1:numel(figCells)
                                        try
                                        close(figCells{figC})
                                        catch
                                        end
                                    end
                                    figCells=[];
                                end
                                PromptHistory.prompt{end+1}=verdict;
                                PromptHistory.composition{end+1}=currentComposition;
                                PromptHistory.curPos=PromptHistory.curPos+1;
                                currentComposition
                                verdict 
                                PromptHistory.curPos
                            end
                            if isempty(verdict)
                                verdict = 0;
                            elseif numel(verdict)>1
                                verdict=real(verdict)+imag(verdict);
                                negativeElementsLocation = find(verdict(:,1)<0);
                                if ~isempty(negativeElementsLocation)
                                    negativeElements = currentComposition(abs(verdict(negativeElementsLocation,1)));
                                    restOfElements = setdiff(currentComposition,negativeElements);
                                    badCompFromInput=[badCompFromInput;allcomb(negativeElements,restOfElements)];
                                    verdict(negativeElementsLocation,:)=[];
                                end
                                if ~isempty(verdict)
                                    badCompFromInput=[badCompFromInput;currentComposition(verdict)];
                                end
                                verdict  = 0;
                            end
                            
                            if isempty(verdict)
                                verdict = 0;
                            end
                            try
                                close(fig_handle);
                            catch
                            end
                            if ~isempty(badCompFromInput)
                                for badC=1:size(badCompFromInput,1)
                                    badcomp = badCompFromInput(badC,:);
                                    testM = squeeze(sum(compositionForest_binArray(:,badcomp),2))==numel(badcomp);
                                    PFArray_purity(testM)=0;
                                end
                            end
                            
                            
                            verdictCount = verdictCount+1;
                        end
                    else
                        verdict = 0;
                    end
                    %                 end
                    clf
                else
                    if any(currentDecisionMatrix(:)==0)
                        verdict = 0;
                    elseif all(currentDecisionMatrix(:)==1)
                        verdict =1;
                    else
                        verdict = GTmode(currentComposition,currentComposition);
                        verdict = prod(verdict(:));
                        verdictCount = verdictCount+1;
                    end
                    
                end
                
                %------------------------------------% plotting finished
                if verdict==1
                    decisionStorage(currentComposition,currentComposition)=1;
                    BFS_purity(node) = 1;
                    PFArray_purity(iterator(node))=1;
                else%-----------%
                    BFS_purity(node) = 0;
                    if ~isempty(badCompFromInput)
                        for badC=1:size(badCompFromInput,1)
                            badcomp = badCompFromInput(badC,:);
                            decisionStorage(badcomp(1),badcomp(2))=0;
                            decisionStorage(badcomp(2),badcomp(1))=0;
                            testM = squeeze(sum(compositionForest_binArray(:,badcomp),2))==numel(badcomp);
                            PFArray_purity(testM)=0;
                        end
                    end
                    if numOfMinusOneEntry==2 && all(currentDecisionMatrix~=0,'all')
                        % if current decison matrix only contains two -1 and (the rest are 1,or it doesn't contains Zero)
                        [M1r,M1c]=find(currentDecisionMatrix==-1); % identify those two -1
                        currentDecisionMatrix(M1r(1),M1c(1))=0; % assign them as zero
                        currentDecisionMatrix(M1r(2),M1c(2))=0;
                        decisionStorage(currentComposition,currentComposition) = currentDecisionMatrix; % update decision Matrix
                        badcomp = currentComposition([M1r(1) M1c(1)]);
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                    end
                    %-----------% verdict ==0 operation finish
                end % verdict processing finished
            else
                if all(currentDecisionMatrix(:)==1)
                    BFS_purity(node) = 1;
                    PFArray_purity(iterator(node))=1;
                else
                    BFS_purity(node) = 0;
                    PFArray_purity(iterator(node))=0;
                end % query finished
            end % previous purity check finished
        end
        
    end
    
    
    
    
    
    UniqueClusterPurity = false(size(UniqueclusterLabel,2),1);
    for u = 1:numel(UniqueClusterPurity)
        UniqueClusterPurity(u)=all(decisionStorage(UniqueclusterLabel(:,u),UniqueclusterLabel(:,u))==1,'all');
        %       if ~UniqueClusterPurity(u)
        %       badcomp = find(UniqueclusterLabel(:,u));
        %       testM = squeeze(sum(compositionForest_binArrayOri(:,badcomp,:),2))==numel(badcomp);
        %       PFArray_purityOri(testM)=0;
        %       end
    end
    
    if sum(~UniqueClusterPurity)>0
        UniqueclusterLabel=UniqueclusterLabel(:, UniqueClusterPurity==0);
        for u = 1:size(UniqueclusterLabel,2)
            if all(min(bsxfun(@minus,UniqueclusterLabel',UniqueclusterLabel(:,u)'),[],2)>0)
                
                badcomp = find(UniqueclusterLabel(:,u));
                testM = squeeze(sum(compositionForest_binArrayOri(:,badcomp,:),2))==numel(badcomp);
                PFArray_purityOri(testM)=0;
            end
        end
    end
    
    
    minExtPurity_partial=double(UniqueClusterPurity(ic));
    minExtPurity = zeros(1,numTree);
    minExtPurity(testD) = minExtPurity_partial;
else
    minExtPurity=zeros(1,numTree);
end