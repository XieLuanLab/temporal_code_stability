function [Z, distCurrentUno] = Link1Dist1onDistanceMatrixClusterLabelMatrix(distMatrix,clu_Label)
numClu = size(clu_Label,2);
distCurrentUno = zeros(numClu,numClu);

for i = 1:numClu-1
                subDistMatrix = distMatrix(clu_Label(:,i),:);

    for j = i+1:numClu
        

            subsubDistMatrix = subDistMatrix(:,clu_Label(:,j));
            allDist = subsubDistMatrix(:);
        
        

    distCurrentUno(i,j)= max(allDist);
    distCurrentUno(j,i)= distCurrentUno(i,j);
    
    end
end

Z = cell(1,2);
for tri = 1:1
   Y = distCurrentUno;
   Y = Y - diag(diag(Y)); % make sure the diagnoal is zero
   yOut = squareform(Y,'tovector');

        Z{1,1} = linkage(yOut,'single');
        Z{1,2} = clu_Label;

end
