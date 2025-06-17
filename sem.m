function final=sem(x)
final=nanstd(x,[],1)./sqrt(sum(~isnan(x),1));
if size(x,1)==1
    'warning'
    final=nanstd(x(:))/sqrt(sum(~isnan(x)));
end
if isempty(x)
final = nan;
end
end