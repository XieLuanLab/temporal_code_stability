function [cluster_output,Nby2] = MutualNearestNeighbour(distMatrix,maskMatrix,sessionNumberList)
distMatrix(distMatrix==0)=nan;% dist must greater than zero, add a small number if that can be zero
newDistMatrix = nan(size(distMatrix));
totalSession = numel(unique(sessionNumberList));
numUnit = size(distMatrix,1);
for thisSes = 1:totalSession-1 % ignoring last two sessions. 
    nextSes = thisSes+1; 
    thisSesCandidateList  = find(sessionNumberList==thisSes ) ;
    nextSesCandidateList  = find(sessionNumberList==nextSes) ;
    newDistMatrix(thisSesCandidateList,nextSesCandidateList)=distMatrix(thisSesCandidateList,nextSesCandidateList); 
end
newDistMatrix(maskMatrix==0)=nan;
Match_Matrix_row = sparse(numUnit,numUnit);
Match_Matrix_col = sparse(numUnit,numUnit);
for i=1:size(Match_Matrix_row,1)
    i;
    if ~isnan(min(newDistMatrix(i,:)))
    [~,LMSE_index]=sort(newDistMatrix(i,:));
  
    selected = LMSE_index(1);
        Match_Matrix_row(i,selected)=1;
    
    end
end

for j=1:size(Match_Matrix_col,2)
    j;
    if ~isnan(min(newDistMatrix(:,j)))
    [~,LMSE_index]=sort(newDistMatrix(:,j));

    selected = LMSE_index(1);
    Match_Matrix_col(selected,j)=1;
    end
end

final_match = Match_Matrix_col.*Match_Matrix_row;
[X,Y]=find(final_match ==1);
Nby2=[X Y];


spreadMatrix  = link_Record(Nby2,totalSession,sessionNumberList);
cluster_output  = -1000*ones(numUnit,1);
for i = 1:size(spreadMatrix,1)
list = spreadMatrix(i,spreadMatrix(i,:)>0);
    cluster_output(list)=i;
end
singletonNum = sum(cluster_output<0);

if singletonNum>0
cluster_output(cluster_output<0) = size(spreadMatrix,1)+(1:singletonNum);
end

