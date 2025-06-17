function [Z, distCurrentTri] = Link4Dist3onDistanceMatrix(distMatrix,clu_Label)
numClu = numel(unique(clu_Label));
distCurrentTri = zeros(numClu,numClu,3);

for i = 1:numClu-1
                subDistMatrix = distMatrix(clu_Label==i,:);

    for j = i+1:numClu
        

            subsubDistMatrix = subDistMatrix(:,clu_Label==j);
            allDist = subsubDistMatrix(:);
        
        
    distCurrentTri(i,j,1)= sum(allDist)/numel(allDist);
    distCurrentTri(j,i,1)=distCurrentTri(i,j,1);
    distCurrentTri(i,j,2)= min(allDist);
    distCurrentTri(j,i,2)=distCurrentTri(i,j,2);
    distCurrentTri(i,j,3)= max(allDist);
    distCurrentTri(j,i,3)=distCurrentTri(i,j,3);
    
    end
end

Z = cell(4,2);
for tri = 1:3
   Y = squeeze(distCurrentTri(:,:,tri));
   Y = Y - diag(diag(Y)); % make sure the diagnoal is zero
   yOut = squareform(Y,'tovector');
% try hierarchical --------------------------------------------------------------------------
switch tri
    case 1
        Z{1,1} = linkage(yOut,'average');
        Z{1,2} = clu_Label;
        Z{4,1} = linkage(yOut,'weighted');
        Z{4,2} = clu_Label;
       
    case 2
        Z{2,1} = linkage(yOut,'single');
        Z{2,2} = clu_Label;

    case 3
        Z{3,1} = linkage(yOut,'complete');
        Z{3,2} = clu_Label;
  
end
end
