function idx=sri(period,nthPeriod,nthIdxWithinPeriod)
% select all elements within the same period from (n-1)th period+1:n the period. 
idx=((nthPeriod-1)*period+1):(nthPeriod*period);
if ~isempty(nthIdxWithinPeriod)
    idx=idx(nthIdxWithinPeriod);
end

