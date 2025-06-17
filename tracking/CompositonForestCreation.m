function compositionForest = CompositonForestCreation(rawZ)
compositionForest = cell(1,size(rawZ,1));

for i =1:numel(compositionForest)

    endNode = size(rawZ{i,1},1)+1 ;
    initialNodeNum =  endNode + size(rawZ{i,1},1) ;

        Z = rawZ{i,1};
        Z(:,4)=( endNode+(1:size(Z,1)))';
        % find parent for each node in Z
        parent = zeros( initialNodeNum, 1);
        for node = 1:initialNodeNum-1
            [r,c] = find(Z(:,1:2)==node);
            parent (node) = Z( r, 4);
        end
        parent (node+1) = 0;
        % initialize tree
        trackingTree = tree(initialNodeNum,0);
        trackingTree.Node = cell(initialNodeNum,1);
        trackingTree.Parent = zeros(initialNodeNum,1);
        
        trackingTree.Node{1} = initialNodeNum;
        trackingTree.Parent(1) = 0;
        
        currentNodes =  trackingTree.Node{1} ;
        currentNodesNewTreeIndex  = 1;
        % construct and finish converting Z variable to a Tree, by recursively
        % adding children to parents.
        NewTreeNodeIndex = 2;
        while sum(ismember(currentNodes,parent))>0 % if there are other Z node that use the currentNodes as parent
            Lia= ismember(currentNodes,parent); % find out currentNode Z index that satisfy this requirement
            validParentNode = currentNodes(Lia); % find out currentNode Z node number that satisty this requirement
            validParentNodeNewTreeIndex  = currentNodesNewTreeIndex(Lia);
            
            currentNodes=[];
            currentNodesNewTreeIndex=[];
            
            for cn =1:numel( validParentNode) % go through each of the valid parent node Z
                nodesUsingValidNodeAsParent = find(parent==validParentNode(cn)); % find all its child nodes
                currentNodes=[currentNodes nodesUsingValidNodeAsParent']; % add child nodes indices Z as new parentNodes
                
                
                for n = nodesUsingValidNodeAsParent'
                    %             trackingTree =  trackingTree.addnode(find(trackingTree==validParentNode(cn)),n);% figure out the index of the parent, add the child node to the parent.
                    
                    trackingTree.Node{NewTreeNodeIndex} = n;
                    trackingTree.Parent(NewTreeNodeIndex) = validParentNodeNewTreeIndex(cn);
                    currentNodesNewTreeIndex = [currentNodesNewTreeIndex NewTreeNodeIndex];
                    NewTreeNodeIndex  = NewTreeNodeIndex+1;
                end
            end
            
        end
        trackingTreeIndiceGroup = trackingTree; % copy the tree
        leafs = trackingTreeIndiceGroup.findleaves;
        leafContent = trackingTreeIndiceGroup.Node(leafs);
        leafContent = cat(1,leafContent{:});
        if size(rawZ,2)>1 % append raw unit content, append parent node as current leaf nodes.
            precluLabel = rawZ{i,2};
            if size(precluLabel,2)==1
                
                
                newNodeStore = cell(numel(unique(precluLabel)),1);
                newParentStore = cell(size(newNodeStore));
                for n = 1:numel(leafs)
                    list = find(precluLabel==leafContent(n));
                    newNodeStore{n}=num2cell(list);
                    newParentStore{n}=ones(numel(list),1)*leafs(n);
                end
                trackingTreeIndiceGroup.Node=[trackingTreeIndiceGroup.Node;cat(1,newNodeStore{:})];
                trackingTreeIndiceGroup.Parent=[trackingTreeIndiceGroup.Parent;cat(1,newParentStore{:})];
            else
                newNodeStore = cell(size(precluLabel,2),1);
                newParentStore = cell(size(newNodeStore));
                for n = 1:numel(leafs)
                    list = find(precluLabel(:,n)); % note this line is CORRECTLY different from its counterpart the previous if 
                    newNodeStore{n}=num2cell(list);
                    newParentStore{n}=ones(numel(list),1)*leafs(n);
                end
                try
                trackingTreeIndiceGroup.Node=[trackingTreeIndiceGroup.Node;cat(2,newNodeStore{:})'];
                catch
                    trackingTreeIndiceGroup.Node=[trackingTreeIndiceGroup.Node;cat(1,newNodeStore{:})];
                end
                trackingTreeIndiceGroup.Parent=[trackingTreeIndiceGroup.Parent;cat(1,newParentStore{:})];

            end
        end
                       joinTree2 = @(x) cat(2,x{:}); % create a function that joins the unit composition of children nodes to form the composition for the parent node.

        compositionForest{i} = trackingTreeIndiceGroup.recursivecumfunCell(@(x) joinTree2(x));  % apply the function recursively
end