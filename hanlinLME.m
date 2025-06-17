function lme=hanlinLME(data,dayStart)
% structure of data : 
[nU,nD]=size(data{1});
nC=numel(data);

DayID_Mat  =  reshape(repmat(repmat([1:nD],nU,1), 1,nC),[],1)-(1-dayStart);
UnitID_Mat  = reshape(repmat(repmat([1:nU]',1,nD),1,nC),[],1);
condID_Mat  = reshape(bsxfun(@plus,zeros(nU*nD,nC),[1:nC]),[],1);
depenent = reshape(cat(2,data{:}),[],1); 

TBL=table(depenent,DayID_Mat,condID_Mat,UnitID_Mat,'VariableNames',{'Y','time','condition','subject'});
TBL.subject=nominal(TBL.subject);TBL.condition=nominal(TBL.condition);

lme1 = fitlme(TBL,['Y~time*condition+(1|subject)']);
lme2 = fitlme(TBL,['Y~time*condition+(1|subject)']);
lme3 = fitlme(TBL,['Y~time*condition+(1|subject)']);
lme4 = fitlme(TBL,['Y~time*condition+(1|subject)']);

lme=[{lme1} {lme2} {lme3} {lme4}];



