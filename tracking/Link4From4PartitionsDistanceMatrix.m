function [Z] = Link4From4PartitionsDistanceMatrix(distMatrix,AllPartition)
Z = cell(4,2);
for tree= 1:4
    clu_Label = AllPartition{tree};
    numClu = numel(unique(clu_Label));
    distCurrentTri = zeros(numClu,numClu);
    
    for i = 1:numClu-1
            subDistMatrix = distMatrix(clu_Label==i,:);

        for j = i+1:numClu
            subsubDistMatrix = subDistMatrix(:,clu_Label==j);
            allDist = subsubDistMatrix(:);
            
            switch tree
                case 1
                    distCurrentTri(i,j)= sum(allDist)/numel(allDist);
                    distCurrentTri(j,i)=distCurrentTri(i,j);
                case 2
                    distCurrentTri(i,j)= min(allDist);
                    distCurrentTri(j,i)=distCurrentTri(i,j);
                case 3
                    distCurrentTri(i,j)= max(allDist);
                    distCurrentTri(j,i)=distCurrentTri(i,j);
                case 4
                    distCurrentTri(i,j)= sum(allDist)/numel(allDist);
                    distCurrentTri(j,i)=distCurrentTri(i,j);    
            end
        end
    end
    Y = distCurrentTri;
    Y = Y - diag(diag(Y)); % make sure the diagnoal is zero
    yOut = squareform(Y,'tovector');
    
    
    
    Z{tree,2} = clu_Label;
    
    switch tree
        case 1
            Z{tree,1} = linkage(yOut,'average');
        case 2
            Z{tree,1} = linkage(yOut,'single');
        case 3
            Z{tree,1} = linkage(yOut,'complete');
        case 4
            Z{tree,1} = linkage(yOut,'weighted');  
    end
    
    

end
