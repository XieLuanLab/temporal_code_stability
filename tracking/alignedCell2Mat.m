function output = alignedCell2Mat(x,value)
IndCount = cellfun(@numel,x);
output = value*ones(max(IndCount),numel(x));
for i = 1:numel(x)
output(1:IndCount(i),i)=x{i};
end
