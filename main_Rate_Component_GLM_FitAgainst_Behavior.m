nVS=3;
% modify stimulus types number (select from 1 ,2 3, 5) to fit for a different type 
% indicated by stimNames, IJO is NI(natrual images)
stimNames = {'DG','SG','RFG','','IJO'};

ustb_fileNameGroups={'pupil_size','limb_speed','pupil_xy1','pupil_xy2','EMG'};
ustb_fileNameGroupsPerShank={'DEL','THE','GAM','PHA'};

codePath=pwd;
basePath=pwd;

addpath(genpath(codePath));

load('nonNoise_whichAnimalLayerShankID7.mat')

trialAverage=false;
timeAverage=true;
ZeroMean=0;
DivideSum=1;
kfold=5;
perDayCVComp=true;
usePCA=false;
daysSelected=1:15;

fit15=true;

kernalWidth=[24];
for ani=[1]
ustart=1;

cd(basePath)

Animal = dir('*visual_FH*');
Animal = Animal([2:4 1 5]);
whichShankIDPerShank={[1:4],[1:3],[1:7],[1:4],[1:7]} % note they are just for counting #elements purpose
TrackSum = cell(1,1);
VSTCB=[];% # of trials ; # of conditions ; # Bins of 0.5ms stim on 
VSTCB.DG= [50 16  1040]; % 520ms
VSTCB.RFG=[30 81  490];  % 245ms
VSTCB.IJO=[60 100 560];  % 280ms
VSTCB.SG= [90 30  520];  % 260ms

cd(Animal(ani).name)
for shank=1:numel(ustb_fileNameGroups)
    ustb_fileName=['USTB_Order_' ustb_fileNameGroups{shank} '.mat']
    shank
    loadingScriptBehavior10Hz
end
TrackSumPerMice=TrackSum;TrackSumPerMice(1)=[];

TrackSum = cell(1,1);
for shank=1:numel(whichShankIDPerShank{ani})
    for feat = 1:numel(ustb_fileNameGroupsPerShank)
    ustb_fileName=['USTB_Order_' ustb_fileNameGroupsPerShank{feat} '_Sh' num2str(shank) '.mat']
    loadingScriptBehavior10Hz
    end
end
TrackSum(1)=[];


%%
temp=[];

VS_Store = cellfun(@(x) eval(['x.' stimNames{nVS}]) ,TrackSumPerMice,'UniformOutput',false);
VS_Store = cat(1,VS_Store {:});
VS_Store = cellfun(@(x) x' ,VS_Store ,'UniformOutput' ,false);
TrackSumPerMice=[];

VS_StorePerSh = cellfun(@(x) eval(['x.' stimNames{nVS}]) ,TrackSum,'UniformOutput',false);
VS_StorePerSh = cat(1,VS_StorePerSh{:});
VS_StorePerSh = cellfun(@(x) x' ,VS_StorePerSh ,'UniformOutput' ,false);

TrackSum=[];

nD=size(VS_Store,2);

nB=size(VS_Store,1);

% special case correct for missing behavioral data for some trials in animal_1, for variable 1,3,4
 loss_trials=isnan(VS_Store{1,9}) | isnan(VS_Store{3,9}) | isnan(VS_Store{4,9});
 temp = (VS_Store{1,8} + VS_Store{1,10})/2 ;
 VS_Store{1,9}(loss_trials)=temp(loss_trials);
 temp = (VS_Store{3,8} + VS_Store{3,10})/2 ;
 VS_Store{3,9}(loss_trials)=temp(loss_trials);
 temp = (VS_Store{4,8} + VS_Store{4,10})/2 ;
 VS_Store{4,9}(loss_trials)=temp(loss_trials);
% ------------------

for feat=[1 2 5]
for d=1:15
    temp=VS_Store{feat,d};
    temp(temp<0.01)=0.01;

switch feat
    case 1
        VS_Store{nB+1,d}=sqrt(temp);
    case 2
        VS_Store{nB+2,d}=log(temp);
    case 5
        VS_Store{nB+3,d}=log(temp);
end

end
end

%--------pupil size sqrt -------
%--------limb speed log -------
%--------- EMG log -------------
%----------pupil segmentation ---------
px=cat(3,VS_Store{3,:});px=px(:); pxz=zscore(px);
py=cat(3,VS_Store{4,:});py=py(:); pyz=zscore(py);

edgeParaList=[0.3:0.01:0.8]
edgeParaStd=zeros(size(edgeParaList));
for e=1:numel(edgeParaList)
edgePara=edgeParaList(e);

pBin=zeros(numel(pxz),5);
pBin(pyz>edgePara & (pyz>pxz & pyz>-pxz),1)=1;
pBin(pxz>edgePara & (pyz<=pxz & pyz>=-pxz) ,2)=1;
pBin(pyz<-edgePara & (pyz<pxz & pyz<-pxz),3)=1;
pBin(pxz<-edgePara & (pyz>=pxz & pyz<=-pxz),4)=1;
pBin(abs(pxz)<=edgePara & abs(pyz)<=edgePara,5)=1;
edgeParaStd(e)=std(sum(pBin));
end
edgePara=edgeParaList(find(edgeParaStd==min(edgeParaStd),1));


pBin=zeros(numel(pxz),1);
pBin(pyz>edgePara & (pyz>pxz & pyz>-pxz))=1;
pBin(pxz>edgePara & (pyz<=pxz & pyz>=-pxz))=2;
pBin(pyz<-edgePara & (pyz<pxz & pyz<-pxz))=3;
pBin(pxz<-edgePara & (pyz>=pxz & pyz<=-pxz))=4;
pBin(abs(pxz)<=edgePara & abs(pyz)<=edgePara)=5;
pBin=reshape(pBin,size(VS_Store{1},1),size(VS_Store{1},2),nD);


nB=size(VS_Store,1);
for d=1:nD
VS_Store{nB+1,d}=pBin(:,:,d);
end
VS_Store(3:4,:)=[]; %remove pupil x and pupil y varible
nB=size(VS_Store,1);

%%
rng(1234)

switch nVS

    case 1
episodeNumWithinBlock=[reshape(repmat([1:25],16,1),[],1);reshape(repmat([1:25],16,1),[],1);];
BlockNum = [ones(400,1);2*ones(400,1)];
totaltrails=800;
nSuper=16;
    case 5
episodeNumWithinBlock=[reshape(repmat([1:15],100,1),[],1);reshape(repmat([1:15],100,1),[],1);reshape(repmat([1:30],100,1),[],1);];
BlockNum = [ones(1500,1);2*ones(1500,1);3*ones(3000,1)];
totaltrails=6000;
nSuper=100;
    case 3
episodeNumWithinBlock=reshape(repmat([1:30],81,1),[],1);
BlockNum = [ones(2430,1);];
totaltrails=2430;
nSuper=81;
    case 2
episodeNumWithinBlock=[reshape(repmat([1:30],30,1),[],1);reshape(repmat([1:30],30,1),[],1);reshape(repmat([1:30],30,1),[],1);];
BlockNum = [ones(900,1);2*ones(900,1);3*ones(900,1)];
totaltrails=2700;
nSuper=30;
end


whichShUnitsOfThisMouse=whichSev7ShankIDPerUnit(whichAnimalPerUnit==ani);


filename=['LDA_' num2str(kernalWidth*2) 'ms_' num2str(kfold) 'pD' num2str(perDayCVComp) 'CV_a_' num2str(ani) 'nVS_' num2str(nVS) '_ZeroMean_' num2str(ZeroMean) '_DivideSum_' num2str(DivideSum) '_PCA' num2str(usePCA) '.mat'];
load([basePath filesep filename],'FR_Store','LDA_15score_Store')


nU=size(FR_Store,1);

LL=zeros(nU,9);
RR=zeros(nU,9);
NZ=cell(nU,9);

if fit15

avail = ~cellfun(@isempty,LDA_15score_Store);% very critical 


availTrain = avail;
availTrain(:,setdiff(1:15,daysSelected)) = false;
for u = ustart:nU
u
if sum(availTrain(u,:))==0
    continue
end
   LFPfeatThisUnit =    VS_StorePerSh(sri(4,whichShUnitsOfThisMouse(u),1:4),:);

for d=1:15
        LFPfeatThisUnit{5,d}=LFPfeatThisUnit{1,d}./LFPfeatThisUnit{2,d};
        LFPfeatThisUnit{6,d}=LFPfeatThisUnit{1,d}./LFPfeatThisUnit{3,d};
        LFPfeatThisUnit{7,d}=LFPfeatThisUnit{2,d}./LFPfeatThisUnit{3,d};
end

temp=FR_Store(u,:);
firingRates=cat(3,temp{availTrain(u,:)}); % remove day axis , as if one day  
Y0 = firingRates(:); %serialization order, all cond from one episode, cond from other episode, cond from other days 

BEH=cell(size(VS_Store,1),1);
for beh=1:size(VS_Store,1)
temp=VS_Store(beh,:);
firingRates=cat(3,temp{availTrain(u,:)}); % remove day axis , as if one day  
BEH{beh} = firingRates(:); %serialization order, all cond from one episode, cond from other episode, cond from other days 
end
BEH{beh} = onehotencode(categorical(BEH{beh}),2); % one hot encode the location zone of pupil. 

for beh=1:size(LFPfeatThisUnit,1)
temp=LFPfeatThisUnit(beh,:);
firingRates=cat(3,temp{availTrain(u,:)}); % remove day axis , as if one day  
BEH{end+1} = firingRates(:); %serialization order, all cond from one episode, cond from other episode, cond from other days 
end

   timeID = repmat(BlockNum,1,sum(availTrain(u,:)));
   BEH{end+1}=timeID(:);
   timeID = repmat(episodeNumWithinBlock,1,sum(availTrain(u,:)));
   BEH{end+1}=timeID(:);
   BEH{end+1}=sqrt(timeID(:));

X0=cat(2,BEH{:});
tic
if var(Y0)==0
    Y0=poissrnd(1,size(Y0));
end
    [B,fitInfo]=lassoglm(X0,Y0,"poisson",'CV',5,'NumLambda',25,'UseCovariance',true);
idxLambda1SE = round((fitInfo.Index1SE+fitInfo.IndexMinDeviance)/2);
mincoefsIncluded = find(B(:,idxLambda1SE));
if isempty(mincoefsIncluded)
    idxLambda1SE=find(any(B>0,1),1,'last');
    mincoefsIncluded = find(B(:,idxLambda1SE));
end
mdl=fitglm(X0(:,mincoefsIncluded),Y0,'Distribution','poisson');
LL(u,8)=fitInfo.Lambda(idxLambda1SE);
RR(u,8)=mdl.Rsquared.Ordinary;
NZ{u,8}=mincoefsIncluded;

% 
for comp=1:3
temp=LDA_15score_Store(u,:);
temp=fillEmptyAsNaN(temp);
temp=cellfun(@(x) x(:,:,comp),temp,'Uni',0);
firingRates=cat(3,temp{availTrain(u,:)}); 
Y0 = firingRates(:); 
if var(Y0)==0
    Y0=randn(size(Y0));
end

[B,fitInfo]=lasso(X0,Y0,'CV',5,'NumLambda',25);
idxLambda1SE = round((fitInfo.Index1SE+fitInfo.IndexMinMSE)/2);
mincoefsIncluded = find(B(:,idxLambda1SE));
if isempty(mincoefsIncluded)
    idxLambda1SE=find(any(B>0,1),1,'last');
    mincoefsIncluded = find(B(:,idxLambda1SE));
end
mdl=fitglm(X0(:,mincoefsIncluded),Y0,'Distribution','normal');
LL(u,comp)=fitInfo.Lambda(idxLambda1SE);
RR(u,comp)=mdl.Rsquared.Ordinary;
NZ{u,comp}=mincoefsIncluded;

end
toc
end

LL15=LL;
RR15=RR;
NZ15=NZ;
end

save([basePath filesep filename],'LL15','RR15','NZ15','-append')


cd ..
end
