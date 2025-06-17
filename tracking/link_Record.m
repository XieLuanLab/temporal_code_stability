function recordMatrix=link_Record(AllMatchRecord,totalSession,lookUpTable)
recordMatrix = NaN(size(AllMatchRecord,1),totalSession); % assume each unit can be matched across all sessions. 
ini = lookUpTable(min(min(AllMatchRecord)),1)-1;
for rec = 1:size(AllMatchRecord,1)
    
   [x,y]=find(recordMatrix== AllMatchRecord(rec,1));
   if isempty(x) % we don't see this record before, write it to the next non-empty line
   [xN,yN]=find(isnan(recordMatrix)~=1);
      if ~isempty(xN)
       recordMatrix(max(xN)+1, lookUpTable(AllMatchRecord(rec,1),1)-ini)=AllMatchRecord(rec,1);
       recordMatrix(max(xN)+1, lookUpTable(AllMatchRecord(rec,2),1)-ini)=AllMatchRecord(rec,2);
      else   
       recordMatrix(0+1, lookUpTable(AllMatchRecord(rec,1),1)-ini)=AllMatchRecord(rec,1);
       recordMatrix(0+1, lookUpTable(AllMatchRecord(rec,2),1)-ini)=AllMatchRecord(rec,2);
       end
   else   % we have seen this record before, write it after last occurance
       recordMatrix(x, lookUpTable(AllMatchRecord(rec,1),1)-ini)=AllMatchRecord(rec,1);
       recordMatrix(x, lookUpTable(AllMatchRecord(rec,2),1)-ini)=AllMatchRecord(rec,2);
   end
end
recordMatrix(isnan(min(recordMatrix')),:)=[];
recordMatrix(isnan(recordMatrix))=-1000;
