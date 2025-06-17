function AllPartition = getPartition(FinalPurityForest,compositionForest)
nTree = numel(compositionForest);
nUnit = max(cell2mat(compositionForest{1}.Node'));



AllPartition  = cell(1,nTree);
for tree= 1:nTree
temp = zeros(1,nUnit);
clu = 1;
    for node = 1:numel(compositionForest{tree}.Node)
   if(FinalPurityForest(node,tree)==true &&  FinalPurityForest(compositionForest{tree}.Parent(node),tree)==false)
   temp(compositionForest{tree}.Node{node})=clu;
   clu = clu +1;
   end
end
AllPartition{tree} = temp';
end


end
