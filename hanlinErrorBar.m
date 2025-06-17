function l=hanlinErrorBar(data,whichColor)
global xaxisshaded
if iscell(data)
   if isempty(xaxisshaded)
    xaxis=1:numel(data);
   else
       xaxis=xaxisshaded;
   end
l=shadedErrorBar(xaxis,cellfun(@nanmean,data),cellfun(@(x) sem(x),data),'lineProps',{'Color',whichColor,'Marker','o','MarkerFaceColor',[1 1 1],'LineWidth',2,'MarkerSize',5});
else 
   if isempty(xaxisshaded)
    xaxis=1:size(data,2);
   else
       xaxis=xaxisshaded;
   end
l=shadedErrorBar(xaxis,nanmean(data,1),sem(data),'lineProps',{'Color',whichColor,'Marker','o','MarkerSize',5,'MarkerFaceColor',[1 1 1],'LineWidth',2});
end
