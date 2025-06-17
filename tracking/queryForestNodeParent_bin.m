
function [largestPureNode,PFArray_purity,decisionStorage,verdictCount] =  queryForestNodeParent_bin(div,compositionForest_binArray,PFArray_purity,PFArray_parent,PFArray_ptr,MinExtNodeIndice,whichTree,decisionStorage,options,GTmode,DataForFeatureAvgPer32days,shankStr)
% This function decide the largest pure node in one tree, starting from
% another node in the same tree, by performing binary search along its path to root

% Input = composition forest, purity forest (parent and purity) , past
% decision, index of node in question, which tree we want search, possible ground truth matrix

% Output = index of the largest pure node, updated purity forest, decision
% matrix, number of feedback required to get decision.

global PromptHistory

verdictCount = 0;
numTree = size(PFArray_purity,2);
compositionTree_binArray = compositionForest_binArray(:,:,whichTree);
purityTree= PFArray_purity(:,whichTree);
purityTreeParent = PFArray_parent(:,whichTree);
largestPureNode = MinExtNodeIndice(whichTree);

endNode = 1 ; % define rootNode as the endNode;
midPointNode = 1;

loop=0
while 1 % enter checking loop
    
    
    %------------------- define break condition
    if (midPointNode==largestPureNode) || (purityTree(largestPureNode)==1 && purityTree(purityTreeParent(largestPureNode))==0)
        % basically when you find the largest pure node, then quit
        % the loop
        largestPureNodeComposition = find(compositionTree_binArray(largestPureNode,:));
        % a verdict 1 at largestPureNode
        % update decision matrix.
        decisionStorage(largestPureNodeComposition,largestPureNodeComposition)=1;
        % update this tree
        
        PFArray_purity(:,whichTree)=purityTree;
        if largestPureNode~=1 % if this largest pure node has a parent, then the parent is bad
            % update other tree
            %-----------% a verdict zero at  parent of largestPureNode.
            
            smallestImpureNodeComposition = find(compositionTree_binArray(purityTreeParent(largestPureNode),:));
            leftChildComp = largestPureNodeComposition;
            rightChildComp = setdiff(smallestImpureNodeComposition ,largestPureNodeComposition);
            leftChildPurity = 1;
            currentDecisionMatrix=decisionStorage(smallestImpureNodeComposition,smallestImpureNodeComposition);
            numOfMinusOneEntry = sum(currentDecisionMatrix(:)==-1);
            crossPurityMatrix = decisionStorage(leftChildComp,rightChildComp);
            crossPurity = all(crossPurityMatrix(:)==1); % whether the two group has the possibility to combine, also need the two group to be pure themselves.
            if all([leftChildPurity==1 crossPurity])
                badcomp = rightChildComp;
                testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                PFArray_purity(testM)=0;
                
                % if cross purity is ok and 1 group is ok, the fault is at the other group.
            else
                badcomp = smallestImpureNodeComposition;
                testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                PFArray_purity(testM)=0;
                
            end
            
            if numOfMinusOneEntry==2 && all(currentDecisionMatrix~=0,'all') % if current decison matrix only contains two -1 and (the rest are 1)
                [M1r,M1c]=find(currentDecisionMatrix==-1); % identify those two -1
                currentDecisionMatrix(M1r(1),M1c(1))=0; % assign them as zero
                currentDecisionMatrix(M1r(2),M1c(2))=0;
                decisionStorage(smallestImpureNodeComposition,smallestImpureNodeComposition) = currentDecisionMatrix; % update decision Matrix
                badcomp = smallestImpureNodeComposition([M1r(1) M1c(1)]);
                testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                PFArray_purity(testM)=0;
                
            end
        end
        
        %-----------% verdict ==0 operation finish
        break;
    end %------------------- define break condition end
    
    
    pathterminate = find(PFArray_ptr(largestPureNode,:)==endNode);
    pathToEnd = PFArray_ptr(largestPureNode,1:pathterminate);
    
    midPointNode = pathToEnd(round(numel(pathToEnd)/div+0.5));
    
    switch purityTree(midPointNode) % binary search,query midPointNode purity
        
        case 0
            if endNode == midPointNode
                midPointNode=largestPureNode;
            end
            endNode = midPointNode;
            continue
        case 1
            largestPureNode = midPointNode;
            continue
        case -1  % prompt user and assign value to purityTree(midPointNode)
            %------------------------------------%obtain info and configure plotting options
            currentNodeChildren = find(purityTreeParent==midPointNode);% --2 identify the childern of the node in the current tre
            leftChild = currentNodeChildren(1);
            try
                rightChild = currentNodeChildren(2); % assume binary tree
            catch
                rightChild=leftChild;
            end
            leftChildComp = sort(find(compositionTree_binArray(leftChild,:)));
            rightChildComp = sort(find(compositionTree_binArray(rightChild,:)));
            if isempty(leftChildComp)
                purityTree(midPointNode)=purityTree(rightChild);
                endNode = rightChild;
                continue
            end
            if isempty(rightChildComp)
                purityTree(midPointNode)=purityTree(leftChild);
                endNode = leftChild;
                continue
            end
            currentComposition = sort(find(compositionTree_binArray(midPointNode,:))); % get unit composition
            curDecisionMatrix = decisionStorage(currentComposition,currentComposition); % obtain subMatrix from decision storage
            try
                if options.requireParentNodeContinuousInTime==1
                    if any(diff(sort(options.dayNumForEachUnit(currentComposition)))>1)
                        curDecisionMatrix(1)=0;
                    end
                end
            catch
            end
            
            
            if isempty(GTmode)
                options.subplotContent = cell(options.subplotRow,options.subplotCol); % specify subplot content
                options.subplotContent{1,1}=leftChildComp;
                options.subplotContent{1,2}=rightChildComp;
                options.subplotContent{2,1}=[leftChildComp(1) rightChildComp(1)];
                options.subplotContent{2,2}=[leftChildComp(end) rightChildComp(end)];
                options.subplotContent{3,1}=[leftChildComp(1) rightChildComp(end)];
                options.subplotContent{3,2}=[leftChildComp(end) rightChildComp(1)];
                options.currentComposition=currentComposition;
                badCompFromInput=[];
                if ~any(curDecisionMatrix(:)==0)
                    
                    
                    if options.defaultPrompt==1
                        t = timer('ExecutionMode', 'singleShot', 'StartDelay', 14+0.04*numel(currentComposition), 'TimerFcn', @pressOneEnter);
                    elseif options.defaultPrompt==2
                        t = timer('ExecutionMode', 'singleShot', 'StartDelay', 14+0.04*numel(currentComposition), 'TimerFcn', @pressTwoEnter);
                    elseif options.defaultPrompt==0
                        t = timer('ExecutionMode', 'singleShot', 'StartDelay', 14+0.04*numel(currentComposition), 'TimerFcn', @pressZeroEnter);
                        
                    end
                    
                    if options.usePromptHistory==1 && PromptHistory.curPos<=numel(PromptHistory.prompt)
                        
                        verdict = PromptHistory.prompt{PromptHistory.curPos};
                        PromptHistory.curPos=PromptHistory.curPos+1;
                        currentComposition
                                verdict 
                                PromptHistory.curPos
                    else
                        for subplotRow = 1:options.subplotRow
                            for subplotCol = 1:options.subplotCol
                                subplot_tight(options.subplotRow,options.subplotCol,(subplotRow-1)*options.subplotCol+subplotCol,[0.01 0.01]);
                                plot_two_units(subplotRow,subplotCol,options,DataForFeatureAvgPer32days)
                                
                            end
                        end
                        displayChoiseButtons()
                        drawnow
                        if sum(options.defaultPrompt==[0 1 2])~=0
                            start(t);
                            verdict = input(['combine P =? 1/yes, 0/no, 2/more, default ' num2str(options.defaultPrompt) ':']);
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
                            verdict = input(['combine P =? 1/yes, 0/no, 2/more, default ' num2str(options.defaultPrompt) ':']);
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
                    
                    
                    if verdict ==2
                        
                        
                        if options.usePromptHistory==1 && PromptHistory.curPos<=numel(PromptHistory.prompt)
                            verdict = PromptHistory.prompt{PromptHistory.curPos};
                            PromptHistory.curPos=PromptHistory.curPos+1;
                            currentComposition
                                verdict 
                                PromptHistory.curPos
                        else
                            uicontrol('Position',[800+950 750-600 60 48],'String','CCG',...
                                'FontSize',14,'Callback', {@CCGgroupPlot,currentComposition,options.dayNumForEachUnit});%'BackgroundColor','r',
                            [fig_handle] = videofig(shankStr,currentComposition,numel(currentComposition), @redrawAcrossDayUnits); % import
                            redrawAcrossDayUnits(1,currentComposition,shankStr);
                            drawnow
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 0.4, 'TimerFcn', @pressEnter);
                            if numel(currentComposition)<=options.minNumFrameToAutoPlay
                                start(t);
                                verdict = input(['combine P =? 1/yes, 0/no :']);
                                if (numel(verdict)==1 && sum(verdict==[0 1 ])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
                                    verdict = input(['past input fault, was wrong enter new one: ']);
                                end
                                stop(t);
                            else
                                verdict = input(['combine P =? 1/yes, 0/no :']);
                                if (numel(verdict)==1 && sum(verdict==[0 1 ])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
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
                        
                        try
                            close(fig_handle);
                        catch
                        end
                    end
                    
                    if ~isempty(badCompFromInput)
                        for badC=1:size(badCompFromInput,1)
                            badcomp = badCompFromInput(badC,:);
                            decisionStorage(badcomp(1),badcomp(2))=0;
                            decisionStorage(badcomp(2),badcomp(1))=0;
                            
                            testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                            PFArray_purity(testM)=0;
                        end
                    end
                    
                    verdictCount = verdictCount+1;
                    %              end
                else
                    verdict = 0;
                end
                clf  %------------------------------------% plotting finished
            else
                CDM = curDecisionMatrix;
                if     any(CDM(:)==0)
                    verdict = 0;
                elseif all(CDM(:)==1)
                    verdict = 1;
                else
                    verdict = GTmode(currentComposition,currentComposition);
                    verdict = prod(verdict(:));
                    verdictCount = verdictCount+1;
                end
                
            end
            purityTree(midPointNode)=verdict;
            %
    end % midPointNodeQueryFinished
end % loop finish function finish