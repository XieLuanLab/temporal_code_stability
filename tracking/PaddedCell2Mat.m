function [outputArg1] = PaddedCell2Mat(x,maxAmount)

outputArg1 = cell2mat(x);
if size(outputArg1,1)<maxAmount
outputArg1 = [outputArg1;repmat(false(1,size(outputArg1,2)),maxAmount-size(outputArg1,1),1)];
end

end

