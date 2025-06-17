function [s1,s2]=takeNpercentSamplesAndTheOtherHalf(temp,trialN,nSampledTrials)
list=randperm(trialN,nSampledTrials);

s1=cellfun(@(x) mean(x(list,:)),temp,'Uni',0 );
s2=cellfun(@(x) mean(x(setdiff(1:trialN,list),:)),temp,'Uni',0 );

end