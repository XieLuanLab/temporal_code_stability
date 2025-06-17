function FinalPurityMatrix   = getFinalPurity(FinalDecisionMatrix,compositionForest)
nTree = numel(compositionForest);
maxNodeNumber = max(cellfun(@(x) numel(x.Parent),compositionForest));
FinalPurityMatrix  = false(maxNodeNumber ,nTree);

for tree=1:nTree
    tree
    thisTreeNodes = compositionForest{tree}.Node;
    nNode = numel(compositionForest{tree}.Node);
for node=1:nNode

        testGroup = thisTreeNodes{node};
        
        FinalPurityMatrix(node,tree)  = all(FinalDecisionMatrix(testGroup,testGroup),'all');
end
end
end
