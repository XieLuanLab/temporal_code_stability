function [designMatrixShuffle,originalList]=shuffleDesignMatrixDecode(designMatrix,shufflePer)
designMatrixShuffle= designMatrix;
for ani=1:size(designMatrix,1)
    numNeuron=size(designMatrix{ani,1},2);
    originalList= 1:numNeuron;
    numShuffle=max(2,ceil(shufflePer(1)/100*numNeuron))*(shufflePer~=0);%!!!!!!!!!!!!
    if numShuffle==0
        return
    end
    for d=1:15
        whichNeuronToShuffle = sort(randperm(numNeuron,numShuffle));
        originalList(whichNeuronToShuffle)=circshift( originalList(whichNeuronToShuffle),1);
        designMatrixShuffle{ani,d}=designMatrixShuffle{ani,d}(:,originalList);
    end
end
