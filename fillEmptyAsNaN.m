function filledCell = fillEmptyAsNaN(inputCell)
y=~cellfun(@isempty,inputCell);
if sum(y(:))>0
firstE = find(y,1);
sizeE = size(inputCell{firstE});
inputCell(~y)={nan(sizeE)};
filledCell=inputCell;
else
    filledCell=inputCell;
end