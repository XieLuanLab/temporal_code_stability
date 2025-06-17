function X=fillmissingHanlin(X,isFR)
%isFR means sample from Poisson
rng(1234)
for row=1:size(X,1)
    attension=isnan( X(row,:) );
    
    if sum(attension)>0
        
        if ~isFR
        X(row,attension)=randn(1,sum(attension))*nanstd(X(row,:))+nanmean(X(row,:));
        else
        X(row,attension)=poissrnd(nanmean(X(row,:)),1,sum(attension));    
        end

    end

end