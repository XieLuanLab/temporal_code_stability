function [MaskTrack_pastdecisions,badCompFromInput]= decodeUserPromptForMaskMatric(numUnit,PromptHistory,options)
% decode prompts and assign maskMtrix accordingly, (later to be incoporated
% into decMatrix(through element wise multiplication) and keep a
% badCompList later to be scanned across PFArray_purity.
MaskTrack_pastdecisions = -1*ones(numUnit,numUnit); % initialize as a N by N -1 array to store current decision.
MaskTrack_pastdecisions = MaskTrack_pastdecisions + eye(size(MaskTrack_pastdecisions))*2; % intialize 1 to diagonal.
badCompFromInput=cell(1);
if options.usePromptHistory ==1 && ~isempty(PromptHistory.prompt{2})
for curPos = 1:PromptHistory.curPos-1
currentComposition=PromptHistory.composition{curPos};
verdict=PromptHistory.prompt{curPos};
contradictoryPromptScanningRange = max(1,curPos-10):min(PromptHistory.curPos-1,curPos+10);
for contraPos = setdiff(contradictoryPromptScanningRange,curPos)
contraComposition = PromptHistory.composition{contraPos};
contraVerdict = PromptHistory.prompt{contraPos};
A=contraComposition;
B=currentComposition;
if isequal(size(A), size(B)) || (isvector(A) && isvector(B) && numel(A) == numel(B))
    if all(A==B,'all')
        if (isequal(contraVerdict,1) && isequal(verdict,1))
        else
    if (isequal(verdict,1) && ~isequal(contraVerdict,2)) || (isequal(contraVerdict,1) && ~isequal(verdict,2)) 
    input(['pleaseCtrlC and change PromptHistory at Pos:' num2str(curPos) ' V ' num2str(contraPos)]);
    end
        end
    end
end
end
if numel(verdict)==1 && sum(verdict==[0 1 2])==0
    verdict = input(['past input fault, was:' num2str(verdict) ' enter new one: ']);
    PromptHistory.prompt{curPos}=verdict;
end
if numel(verdict)==1 && verdict==2 % skip ask for more information prompt 
continue
end
if numel(verdict)==1 && verdict==1 % user confirmed merge 
    MaskTrack_pastdecisions(currentComposition,currentComposition)=1;
continue
end

if isempty(verdict) || (numel(verdict)==1 && verdict==0) % user rejected merge
    if numel(currentComposition)==2
    MaskTrack_pastdecisions(currentComposition(1),currentComposition(2))=0;
    MaskTrack_pastdecisions(currentComposition(2),currentComposition(1))=0;
    end
        badCompFromInput{end+1}=currentComposition;

continue
end
[ri,ci]=find(imag(verdict)~=0); % user has uncertain merge, remove them from consideration 
if ~isempty(ri)
   verdict(ri,:)=[];
   if isempty(verdict)
       continue
   end
   
end

    negativeElementsLocation = find(verdict(:,1)<0); % true negative "should not merge" pairs. 
    if ~isempty(negativeElementsLocation)
        negativeElements = currentComposition(abs(verdict(negativeElementsLocation,1)));
        restOfElements = setdiff(currentComposition,negativeElements);
        BCFI=allcomb(negativeElements,restOfElements);
        for badE=1:size(BCFI,1)
                MaskTrack_pastdecisions(BCFI(badE,1),BCFI(badE,2))=0;
                MaskTrack_pastdecisions(BCFI(badE,2),BCFI(badE,1))=0;
        end
        verdict(negativeElementsLocation,:)=[];
    end
    if ~isempty(verdict)
        BCFI=currentComposition(verdict);
        for badE=1:size(BCFI,1)
                MaskTrack_pastdecisions(BCFI(badE,1),BCFI(badE,2))=0;
                MaskTrack_pastdecisions(BCFI(badE,2),BCFI(badE,1))=0;
        end
    end
end
end

