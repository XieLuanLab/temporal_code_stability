
function [minExtPurity,PFArray_purity,decisionStorage,verdictCount,BadCompInputStore] = queryForestNodeChildren_bin(compositionForest_binArray,PFArray_purity,PFArray_parent,MinExtNodeIndicesPerTree,decisionStorage,options,GTmode,DataForFeatureAvgPer32days,shankStr)
% This function decide the purity of one node for each tree inside the
% whole forest.

% Input = composition forest, purity forest (parent and purity) , past
% decision, indices of node in question, possible ground truth matrix

% Output = purity of the nodes in question, updated purity forest, decision
% matrix, number of feedback required to get decision.

% This is a core process in Alg1 in the Zhu, H., Li, X., Sun, L., He, F., Zhao, Z., Luan, L., ... & Xie, C. (2020). 
% Clustering with fast, automated and reproducible assessment applied to longitudinal neural tracking. arXiv preprint arXiv:2003.08533.

% This code also examplify the quote "Purity satisfies important inequalities
% imposed by set inclusion" by simultaneously
% eliminating "bad compositions" across all trees inside the forest. This
% is a key innovation for multi-tree setup. search "badcomp" for more
% information
global PromptHistory

BadCompInputStore = [];
verdictCount = 0;
numTree = size(PFArray_purity,2);
minExtPurity = -1*ones(1,numTree); % initilize
restart=0;
while any(minExtPurity==-1)
    restart=restart+1;
    if restart==3
        break
    end
    for tr = 1:numTree % --1 loop through each tree
        
        minExtPurity(tr) = PFArray_purity(MinExtNodeIndicesPerTree(tr),tr);% are there previous decisions available.
        badCompFromInput=[];
        if PFArray_purity(MinExtNodeIndicesPerTree(tr),tr)==-1 % when previous label doesn't exist.
            
            currentComposition = find(compositionForest_binArray(MinExtNodeIndicesPerTree(tr),:,tr)); % get unit composition
            %         if all(ismember([22 105],currentComposition)) || all(ismember([22 224],currentComposition))
            %             'interesting'
            %         end
            
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
                currentNodeChildren = find(PFArray_parent(:,tr)==MinExtNodeIndicesPerTree(tr));% --2 identify the childern of the node in the current tre
                leftChild = currentNodeChildren(1);
                rightChild = currentNodeChildren(2); % assume binary tree
                leftChildComp = sort(find(compositionForest_binArray(leftChild,:,tr)));
                rightChildComp = sort(find(compositionForest_binArray(rightChild,:,tr)));
                
                DL = decisionStorage(leftChildComp,leftChildComp);
                DL=DL(:);
                if all(DL==1)
                    leftChildPurity=1;
                elseif any(DL==0)
                    leftChildPurity=0;
                elseif numel(currentNodeChildren)>2
                    leftChildPurity=-1;
                else
                    leftChildPurity=PFArray_purity(leftChild,tr);
                end
                
                DR = decisionStorage(rightChildComp,rightChildComp);
                DR=DR(:);
                if all(DR==1)
                    rightChildPurity=1;
                elseif any(DR==0)
                    rightChildPurity=0;
                elseif numel(currentNodeChildren)>2
                    rightChildPurity=-1;
                else
                    rightChildPurity=PFArray_purity(rightChild,tr);
                end
                
                
                if numel(currentNodeChildren)>2
                    leftChild=MinExtNodeIndicesPerTree(tr);
                    rightChild= MinExtNodeIndicesPerTree(tr);
                    leftChildComp = currentComposition;
                    leftChildPurity =-1;
                    rightChildComp = currentComposition;
                    rightChildPurity =-1;
                end
                
                crossPurityMatrix = decisionStorage(leftChildComp,rightChildComp);
                crossPurity = all(crossPurityMatrix(:)==1); % whether the two group has the possibility to combine, also need the two group to be pure themselves.
                % might be better to assign -1
                
                
                if isempty(GTmode) % obtain feedback from users either by displaying two units 's waveform or by looking at ground truth
                    options.currentComposition=currentComposition;
                    options.subplotContent{1,1}=leftChildComp;
                    options.subplotContent{1,2}=rightChildComp;
                    options.subplotContent{2,1}=[leftChildComp(1) rightChildComp(1)];
                    options.subplotContent{2,2}=[leftChildComp(end) rightChildComp(end)];
                    options.subplotContent{3,1}=[leftChildComp(1) rightChildComp(end)];
                    options.subplotContent{3,2}=[leftChildComp(end) rightChildComp(1)];
                    if ~any(currentDecisionMatrix(:)==0)
                        
                        if options.defaultPrompt==1
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 14+numel(currentComposition)*0.04, 'TimerFcn', @pressOneEnter);
                        elseif options.defaultPrompt==2
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 1.5, 'TimerFcn', @pressTwoEnter);
                        elseif options.defaultPrompt==0
                            t = timer('ExecutionMode', 'singleShot', 'StartDelay', 100, 'TimerFcn', @pressZeroEnter);
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
                                    verdict = input(['past input fault, was wrong  enter new one: ']);
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
                            %                     plot_perDay_units(options);
                            
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
                                    verdict = input(['combine=? 1/yes, 0/no :']);
                                    if (numel(verdict)==1 && sum(verdict==[0 1])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
                                        verdict = input(['past input fault, was wrong  enter new one: ']);
                                    end
                                    stop(t);
                                else
                                    verdict = input(['combine=? 1/yes, 0/no :']);
                                    if (numel(verdict)==1 && sum(verdict==[0 1])==0) || (numel(verdict)>1 && max(real(verdict)+imag(verdict),[],'all')>numel(currentComposition))
                                        verdict = input(['past input fault, was wrong  enter new one: ']);
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
                                    testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                                    PFArray_purity(testM)=0;
                                end
                            end
                            
                            
                            verdictCount = verdictCount+1;
                        end
                        DL = decisionStorage(leftChildComp,leftChildComp);
                        DL=DL(:);
                        if all(DL==1)
                            leftChildPurity=1;
                        elseif any(DL==0)
                            leftChildPurity=0;
                        elseif numel(currentNodeChildren)>2
                            leftChildPurity=-1;
                        else
                            leftChildPurity=PFArray_purity(leftChild,tr);
                        end
                        
                        DR = decisionStorage(rightChildComp,rightChildComp);
                        DR=DR(:);
                        if all(DR==1)
                            rightChildPurity=1;
                        elseif any(DR==0)
                            rightChildPurity=0;
                        elseif numel(currentNodeChildren)>2
                            rightChildPurity=-1;
                        else
                            rightChildPurity=PFArray_purity(rightChild,tr);
                        end
                        if (verdict ==0 && numel(currentNodeChildren)==2) && leftChildPurity==-1
                            if numel(leftChildComp)>1
                                
                                if options.usePromptHistory==1 && PromptHistory.curPos<=numel(PromptHistory.prompt)
                                    leftChildPurity = PromptHistory.prompt{PromptHistory.curPos};
                                    PromptHistory.curPos=PromptHistory.curPos+1;
                                    leftChildComp
                                leftChildPurity 
                                PromptHistory.curPos
                                else
                                    uicontrol('Position',[800+950 750-600 60 48],'String','CCG',...
                                        'FontSize',14,'Callback', {@CCGgroupPlot,leftChildComp,options.dayNumForEachUnit});%'BackgroundColor','r',
                                    [fig_handle] = videofig(shankStr,leftChildComp,numel(leftChildComp), @redrawAcrossDayUnits); % import
                                    redrawAcrossDayUnits(1,leftChildComp,shankStr);
                                    drawnow
                                    t = timer('ExecutionMode', 'singleShot', 'StartDelay', 0.4, 'TimerFcn', @pressEnter);
                                    if numel(leftChildComp)<=options.minNumFrameToAutoPlay
                                        start(t);
                                        leftChildPurity = input(['left=? 1/yes, 0/no :']);
                                        if numel(leftChildPurity)==1 && sum(leftChildPurity==[0 1])==0 || (numel(leftChildPurity)>1 && max(real(leftChildPurity)+imag(leftChildPurity),[],'all')>numel(leftChildComp))
                                           leftChildPurity = input(['past input fault, was wrong enter new one: ']);
                                        end
                                        stop(t);
                                    else
                                        leftChildPurity = input(['left=? 1/yes, 0/no :']);
                                        if numel(leftChildPurity)==1 && sum(leftChildPurity==[0 1])==0 || (numel(leftChildPurity)>1 && max(real(leftChildPurity)+imag(leftChildPurity),[],'all')>numel(leftChildComp))
                                            leftChildPurity = input(['past input fault, was wrong enter new one: ']);
                                        end
                                        
                                    end
                                
                                    PromptHistory.prompt{end+1}=leftChildPurity;
                                    PromptHistory.composition{end+1}=leftChildComp;
                                    PromptHistory.curPos=PromptHistory.curPos+1;
                                    leftChildComp
                                leftChildPurity
                                PromptHistory.curPos
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
                                end
                                if isempty(leftChildPurity)
                                    leftChildPurity = 0;
                                elseif numel(leftChildPurity)>1
                                    leftChildPurity=real(leftChildPurity)+imag(leftChildPurity);
                                    negativeElementsLocation = find(leftChildPurity(:,1)<0);
                                    if ~isempty(negativeElementsLocation)
                                        negativeElements = leftChildComp(abs(leftChildPurity(negativeElementsLocation,1)));
                                        restOfElements = setdiff(leftChildComp,negativeElements);
                                        badCompFromInput=[badCompFromInput;allcomb(negativeElements,restOfElements)];
                                        leftChildPurity(negativeElementsLocation,:)=[];
                                    end
                                    if ~isempty(leftChildPurity)
                                        badCompFromInput=[badCompFromInput;leftChildComp(leftChildPurity)];
                                    end
                                    leftChildPurity = 0;
                                end
                                
                                
                                try
                                    close(fig_handle);
                                catch
                                end
                                verdictCount = verdictCount+1;
                            else
                                leftChildPurity=1;
                                
                            end
                        end
                        
                        if (verdict ==0 && numel(currentNodeChildren)==2) && rightChildPurity==-1
                            if numel(rightChildComp)>1
                                
                                if options.usePromptHistory==1 && PromptHistory.curPos<=numel(PromptHistory.prompt)
                                    rightChildPurity = PromptHistory.prompt{PromptHistory.curPos};
                                    PromptHistory.curPos=PromptHistory.curPos+1;
                                    rightChildComp
                                    rightChildPurity
                                     PromptHistory.curPos
                                else
                                    uicontrol('Position',[800+950 750-600 60 48],'String','CCG',...
                                        'FontSize',14,'Callback', {@CCGgroupPlot,rightChildComp,options.dayNumForEachUnit});%'BackgroundColor','r',
                                    
                                    [fig_handle] = videofig(shankStr,rightChildComp,numel(rightChildComp), @redrawAcrossDayUnits); % import
                                    redrawAcrossDayUnits(1,rightChildComp,shankStr);
                                    drawnow
                                    t = timer('ExecutionMode', 'singleShot', 'StartDelay', 0.4, 'TimerFcn', @pressEnter);
                                    if numel(rightChildComp)<=options.minNumFrameToAutoPlay
                                        start(t);
                                        rightChildPurity = input(['right=? 1/yes, 0/no :']);
                                        if numel(rightChildPurity)==1 && sum(rightChildPurity==[0 1])==0 || (numel(rightChildPurity)>1 && max(real(rightChildPurity)+imag(rightChildPurity),[],'all')>numel(rightChildComp))
                                            rightChildPurity = input(['past input fault, was wrong enter new one: ']);
                                        end
                                        stop(t);
                                    else
                                        rightChildPurity = input(['right=? 1/yes, 0/no :']);
                                        if numel(rightChildPurity)==1 && sum(rightChildPurity==[0 1 ])==0  || (numel(rightChildPurity)>1 && max(real(rightChildPurity)+imag(rightChildPurity),[],'all')>numel(rightChildComp))
                                            rightChildPurity = input(['past input fault, was wrong enter new one: ']);
                                        end
                                        
                                    end
                                   
                                    PromptHistory.prompt{end+1}=rightChildPurity;
                                    PromptHistory.composition{end+1}=rightChildComp;
                                    PromptHistory.curPos=PromptHistory.curPos+1;
                                    rightChildComp
                                rightChildPurity
                                PromptHistory.curPos
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
                                end
                                if isempty(rightChildPurity)
                                    rightChildPurity = 0;
                                elseif numel(rightChildPurity)>1
                                    rightChildPurity=real(rightChildPurity)+imag(rightChildPurity);
                                    negativeElementsLocation = find(rightChildPurity(:,1)<0);
                                    if ~isempty(negativeElementsLocation)
                                        negativeElements = rightChildComp(abs(rightChildPurity(negativeElementsLocation,1)));
                                        restOfElements = setdiff(rightChildComp,negativeElements);
                                        badCompFromInput=[badCompFromInput;allcomb(negativeElements,restOfElements)];
                                        rightChildPurity(negativeElementsLocation,:)=[];
                                    end
                                    if ~isempty(rightChildPurity)
                                        badCompFromInput=[badCompFromInput;rightChildComp(rightChildPurity)];
                                    end
                                    
                                    rightChildPurity = 0;
                                end
                                
                                
                                try
                                    close(fig_handle);
                                catch
                                end
                                verdictCount = verdictCount+1;
                            else
                                rightChildPurity=1;
                            end
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
                if leftChildPurity==1
                    decisionStorage(leftChildComp,leftChildComp)=1;
                    PFArray_purity(leftChild,tr)=1;
                end
                if rightChildPurity==1
                    decisionStorage(rightChildComp,rightChildComp)=1;
                    PFArray_purity(rightChild,tr)=1;
                end
                if verdict==1
                    decisionStorage(currentComposition,currentComposition)=1;
                    minExtPurity(tr) = 1;
                    PFArray_purity(MinExtNodeIndicesPerTree(tr),tr)=1;
                elseif verdict == 5
                    continue
                else%-----------%
                    minExtPurity(tr) = 0;
                    
                    if ~isempty(badCompFromInput)
                        for badC=1:size(badCompFromInput,1)
                            badcomp = badCompFromInput(badC,:);
                            decisionStorage(badcomp(1),badcomp(2))=0;
                            decisionStorage(badcomp(2),badcomp(1))=0;
                            testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                            PFArray_purity(testM)=0;
                        end
                    end
                    
                    if all([leftChildPurity==1 crossPurity==1]) || all([rightChildPurity==1 crossPurity==1])
                        badcomp = setdiff(currentComposition,[leftChildComp*(leftChildPurity==1) rightChildComp*(rightChildPurity==1)]);
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                    elseif leftChildPurity==0 && rightChildPurity==1 % left Child is not ok
                        badcomp = leftChildComp;
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                    elseif leftChildPurity==1 && rightChildPurity==0 % right Child is not ok
                        badcomp = rightChildComp;
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                    elseif leftChildPurity==0 && rightChildPurity==0 % No Children is ok
                        badcomp = rightChildComp;
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                        badcomp = leftChildComp;
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                    else
                        badcomp = currentComposition;
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                        
                    end
                    
                    
                    
                    if numOfMinusOneEntry==2 && all(currentDecisionMatrix~=0,'all')
                        % if current decison matrix only contains two -1 and (the rest are 1,or it doesn't contains Zero)
                        [M1r,M1c]=find(currentDecisionMatrix==-1); % identify those two -1
                        currentDecisionMatrix(M1r(1),M1c(1))=0; % assign them as zero
                        currentDecisionMatrix(M1r(2),M1c(2))=0;
                        decisionStorage(currentComposition,currentComposition) = currentDecisionMatrix; % update decision Matrix
                        badcomp = currentComposition([M1r(1) M1c(1)]);
                        testM = squeeze(sum(compositionForest_binArray(:,badcomp,:),2))==numel(badcomp);
                        PFArray_purity(testM)=0;
                    end
                    
                    
                    %-----------% verdict ==0 operation finish
                    
                end % verdict processing finished
            else
                if all(currentDecisionMatrix(:)==1)
                    minExtPurity(tr) = 1;
                    PFArray_purity(MinExtNodeIndicesPerTree(tr),tr)=1;
                else
                    minExtPurity(tr) = 0;
                    PFArray_purity(MinExtNodeIndicesPerTree(tr),tr)=0;
                    
                end % query finished
            end % previous purity check finished
        end % current tree process finished
        BadCompInputStore=[BadCompInputStore;badCompFromInput];
    end % function finish
end  % no uncertainty remains