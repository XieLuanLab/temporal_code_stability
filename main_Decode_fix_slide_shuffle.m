%% 
basePath = pwd;
addpath(genpath(basePath));

for nVS=[1 2 3 5]

daysSelected=1:7;
addpath(genpath(pwd));
load('nonNoise_whichAnimalLayerShankID7.mat')
stimNames = {'DG','SG','RFG','','IJO'};
cd([basePath filesep stimNames{nVS}])
load('VS_Store2025_selectedOnlyRefit.mat')
load([basePath filesep stimNames{nVS} filesep 'pFinalFinalSelect.mat'])

VS_Store(avail==0)={nan(size(VS_Store{1}))};
VS_Store=cellfun(@(x) x',VS_Store,'UniformOutput',0);


load('umapWorstDistance.mat')
umapWorstDistance=nan(size(avail,1),1);
umapWorstDistance(SU2_acrossSes & sum(avail,2)>=3)=cat(1,maxInterDayDistStore{:});
umapThres=85;

availMoreThan = sum(avail,2)>=15 & all(avail(:,(max(daysSelected)+1):end),2);

if daysSelected(1)>2
    availMoreThan = availMoreThan & all(avail(:,daysSelected),2);
end

finalSelect= SU2_acrossSes  & availMoreThan;

rng(1234)
switch nVS
    case 1
groundTruthStimPattern_linear=reshape(repmat([1:16],50,1)',[],1);
nSuper=16;
totaltrails=800;
    case 5
groundTruthStimPattern_linear=reshape(repmat([1:100],60,1)',[],1);
nSuper=100;
totaltrails=6000;
case 3
groundTruthStimPattern_linear=reshape(repmat([1:81],30,1)',[],1);
nSuper=81;
totaltrails=2430;
    case 2
groundTruthStimPattern_linear=reshape(repmat([1:30],90,1)',[],1);
nSuper=30 ;     
totaltrails=2700;
end
dayID = repmat(daysSelected,totaltrails,1);
dayID=dayID(:);
train_trials=arrayfun(@(x) sri(totaltrails,x,1:floor(totaltrails*(numel(daysSelected)-1)/numel(daysSelected))),1:numel(daysSelected),'UniformOutput',false);
train_trials=cat(2,train_trials{:});
valid_trails=setdiff(1:(totaltrails*numel(daysSelected)),train_trials);


groundTruthStimPattern_linear_test=groundTruthStimPattern_linear;
gTSPLA=1:numel(groundTruthStimPattern_linear_test);

ErrorArrayAnimal=cell(1,5);
decoderInfoAnimal=cell(3,5);
a=7.4;
[x,y]=meshgrid(-4:4,-4:4);
x=x*3.7;
y=flip(y)*3.7;
d=20.5;
azim = 45-atand(x/d);
elev = atand(y./sqrt(x.^2+d^2));
azimSerial = azim';azimSerial =azimSerial(:);
elevSerial= elev'; elevSerial=elevSerial(:);
AzEv=[azimSerial elevSerial];
AngleErrorMatrix=zeros(81,81);

for i =1:81
    for j=1:81
        if j~=i
    AzEl1 = AzEv(i,:);
AzEl2 = AzEv(j,:);
a = [[cosd(AzEl1(1)) sind(AzEl1(1))] * cosd(AzEl1(2)) sind(AzEl1(2))];
b = [[cosd(AzEl2(1)) sind(AzEl2(1))] * cosd(AzEl2(2)) sind(AzEl2(2))];
AngleErrorMatrix(i,j) = acosd(dot(a,b));
        end
        end
end

%% 
for ani=1:5
selUid=find(finalSelect & whichAnimalPerUnit==ani & umapWorstDistance<prctile(umapWorstDistance(whichAnimalPerUnit==ani),umapThres));

decoderInfoAnimal{1,ani}=selUid;
fc=VS_Store(selUid,:);

X=[];
Xnew=[];
for dayS=1:numel(daysSelected)
    X=[X cell2mat(cellfun(@(x) reshape(x,[],1)',fc(:,daysSelected(dayS)),'UniformOutput', false))];
end
YY = repmat(groundTruthStimPattern_linear,numel(daysSelected),1);


 filledX=fillmissingHanlin(X,true);
 badrowsX=var(filledX,[],2)<0.001;
 X=X(~badrowsX,:);

XXX=X(:,valid_trails);YYY=YY(valid_trails);X=X(:,train_trials);YY=YY(train_trials);

data_X=fillmissingHanlin(X,true)';
       [Mdl] = fitcdiscr(data_X,YY,...
    'OptimizeHyperparameters',{'Gamma'},...
    'HyperparameterOptimizationOptions',struct('ShowPlots',false, ...
    'UseParallel',false,'KFold',5,...
    'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30,'Repartition',true));   
    Mdl.Gamma=Mdl.HyperparameterOptimizationResults.XAtMinObjective.Gamma;

decoderInfoAnimal{2,ani}=Mdl.Gamma;



for day =7:15 % 
    day
    Xnew=cell2mat(cellfun(@(x) reshape(x,[],1)', fc(:,day),'UniformOutput', false)); % 
    if day==7
        Xnew=XXX;
    end

    Xnew=Xnew(~badrowsX,:);
    badrows=sum(isnan(Xnew),2)<0; % no row is bad
    Xnew=Xnew(~badrows,:);


data_X=fillmissingHanlin(X(~badrows,:),true)';

Xnew=fillmissingHanlin([Xnew X(~badrows,:)],true);

if day==7
    Xnew=Xnew(:,1:numel(YYY));
else
    Xnew=Xnew(:,1:numel(groundTruthStimPattern_linear));
end

    

    Mdl = fitcdiscr(data_X,YY,'Gamma',Mdl.Gamma);

    [maxI,IdxProb]=predict(Mdl,Xnew');maxI=maxI(:);maxI=maxI';

  

    decodedAngle = (360*(maxI-1)/nSuper)';
    ydegreeNew= (groundTruthStimPattern_linear_test-1)*22.5; % ground truth new data
      if day==7
                        ydegreeNew=(YYY-1)*22.5;
                        end
    error = abs(decodedAngle - ydegreeNew);
    error(error>180)=360-error(error>180);
    ErrorArray{1,day}=error;
    
    if nVS==1
    angleOneHot=[0:22.5:337.5]';
    angleOneHot=cosd(angleOneHot)+1j*sind(angleOneHot);
    weightedAngle=IdxProb*angleOneHot;
    decodedAngleW = wrapTo360(rad2deg(angle(weightedAngle)));
    error = abs(decodedAngleW - ydegreeNew);
    error(error>180)=360-error(error>180);
    ErrorArray{3,day}=error;
    end


    if nVS==3 && day==7 
        error=AngleErrorMatrix(sub2ind([size(AngleErrorMatrix,1), size(AngleErrorMatrix,2)],maxI',YYY));
        ErrorArray{1,day}=error;
    elseif nVS==3 && day~=7
        error=AngleErrorMatrix(sub2ind([size(AngleErrorMatrix,1),size(AngleErrorMatrix,2)],maxI',groundTruthStimPattern_linear_test));
        ErrorArray{1,day}=error;
    end


      if day==7
                        ErrorArray{2,day}=maxI'==YYY;
      else
          ErrorArray{2,day}=maxI'==groundTruthStimPattern_linear_test;
      end


    if day~=7
    [mean(error) mean(maxI'==groundTruthStimPattern_linear_test)]
    else
        [mean(error) mean(maxI'==YYY)]
    end

end
ErrorArrayAnimal{ani}=ErrorArray;
end
save(['ErrorArrayAnimal_ratesOnly'  '.mat'],'ErrorArrayAnimal','decoderInfoAnimal')
%% shuffle 
addpath(genpath(basePath));

shufflePer=10;

load('nonNoise_whichAnimalLayerShankID7.mat')
 stimNames = {'DG','SG','RFG','','IJO'};
cd([basePath filesep stimNames{nVS}])
load('VS_Store2025_selectedOnlyRefit.mat')
VS_Store(avail==0)={nan(size(VS_Store{1}))};
VS_Store=cellfun(@(x) x',VS_Store,'UniformOutput',0);
load('umapWorstDistance.mat')
umapWorstDistance=nan(size(avail,1),1);
umapWorstDistance(SU2_acrossSes & sum(avail,2)>=3)=cat(1,maxInterDayDistStore{:});
umapThres=85;

availMoreThan = sum(avail,2)>=15 & sum(avail(:,end-2:end),2)>0;
finalSelect= SU2_acrossSes  & availMoreThan;% 
daysSelected=1:7;

rng(1234)
switch nVS
    case 1
groundTruthStimPattern_linear=reshape(repmat([1:16],50,1)',[],1);
nSuper=16
totaltrails=800
    case 5
groundTruthStimPattern_linear=reshape(repmat([1:100],60,1)',[],1);
nSuper=100
totaltrails=6000
case 3
groundTruthStimPattern_linear=reshape(repmat([1:81],30,1)',[],1);
nSuper=81
totaltrails=2430
    case 2
groundTruthStimPattern_linear=reshape(repmat([1:30],90,1)',[],1);
nSuper=30      
totaltrails=2700
end
train_trials=arrayfun(@(x) sri(totaltrails,x,1:floor(totaltrails*6/7)),1:7,'UniformOutput',false);
train_trials=cat(2,train_trials{:});

valid_trails=setdiff(1:(totaltrails*7),train_trials);

dayID = repmat(daysSelected,totaltrails,1);
dayID=dayID(:);




ErrorArrayAnimal=cell(1,5);
decoderInfoAnimal=cell(3,5);

for ani=1:5

    selUid=find(finalSelect & whichAnimalPerUnit==ani & umapWorstDistance<prctile(umapWorstDistance(whichAnimalPerUnit==ani),umapThres));
    decoderInfoAnimal{1,ani}=selUid;
    fc=VS_Store(selUid,:);


    designMatrix=cell(1,15);
    for aaa=1:1
        for ddd=1:15
            designMatrix{aaa,ddd}=fc(:,ddd)';
        end
    end
    [designMatrixShuffle,originalList]=shuffleDesignMatrixDecode(designMatrix,shufflePer);
    designMatrixShuffle=cellfun(@(x) x',designMatrixShuffle,'UniformOutput',false);
    designMatrixShuffle=cat(1,designMatrixShuffle{:});
    designMatrixShuffle=reshape(designMatrixShuffle,[],15);
    fc=designMatrixShuffle;designMatrixShuffle=[];designMatrix=[];


    X=[];
    Xnew=[];
    for dayS=1:numel(daysSelected)
        X=[X cell2mat(cellfun(@(x) reshape(x,[],1)',fc(:,daysSelected(dayS)),'UniformOutput', false))];
    end
    YY = repmat(groundTruthStimPattern_linear,numel(daysSelected),1);
    XXX=X(:,valid_trails);YYY=YY(valid_trails);X=X(:,train_trials);YY=YY(train_trials);


       [Mdl] = fitcdiscr(fillmissingHanlin(X)',YY,...
    'OptimizeHyperparameters',{'Gamma'},...
    'HyperparameterOptimizationOptions',struct('ShowPlots',false, ...
    'UseParallel',false,'KFold',5,...
    'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30,'Repartition',true));   
     Mdl.Gamma=Mdl.HyperparameterOptimizationResults.XAtMinObjective.Gamma;

decoderInfoAnimal{2,ani}=Mdl.Gamma;

for day =7:15 % 
    day
    Xnew=cell2mat(cellfun(@(x) reshape(x,[],1)', fc(:,day),'UniformOutput', false)); % 
    if day==7
        Xnew=XXX;
    end
    badrows=sum(isnan(Xnew),2)>0;
    Xnew=Xnew(~badrows,:);

    Mdl = fitcdiscr(fillmissingHanlin(X(~badrows,:))',YY,'Gamma',Mdl.Gamma);

    Xnew=fillmissingHanlin([Xnew X(~badrows,:)]);

if day==7
    Xnew=Xnew(:,1:numel(YYY));
else
    Xnew=Xnew(:,1:numel(groundTruthStimPattern_linear));
end
    maxI=predict(Mdl,Xnew')';
    decodedAngle = (360*(maxI-1)/nSuper)';
    ydegreeNew= (groundTruthStimPattern_linear-1)*22.5; % ground truth new data
      if day==7
                        ydegreeNew=(YYY-1)*22.5;
                        end
    error = abs(decodedAngle - ydegreeNew);
    error(error>180)=360-error(error>180);
    ErrorArray{1,day}=error;

    
    
      if day==7
                        ErrorArray{2,day}=maxI'==YYY;
      else
          ErrorArray{2,day}=maxI'==groundTruthStimPattern_linear;
      end

    if day~=7
    [mean(error) mean(maxI'==groundTruthStimPattern_linear)]
    else
        [mean(error) mean(maxI'==YYY)]
    end
end
ErrorArrayAnimal{ani}=ErrorArray;

end
save(['ErrorArrayAnimal_ratesOnly_Shuf'  '.mat'],'ErrorArrayAnimal','decoderInfoAnimal')
%% Sliding 
addpath(genpath(basePath));
load('nonNoise_whichAnimalLayerShankID7.mat')
 stimNames = {'DG','SG','RFG','','IJO'};
cd([basePath filesep stimNames{nVS}])
load('VS_Store2025_selectedOnlyRefit.mat')
VS_Store(avail==0)={nan(size(VS_Store{1}))};
VS_Store=cellfun(@(x) x',VS_Store,'UniformOutput',0);

load('umapWorstDistance.mat')
umapWorstDistance=nan(size(avail,1),1);
umapWorstDistance(SU2_acrossSes & sum(avail,2)>=3)=cat(1,maxInterDayDistStore{:});
umapThres=85;

availMoreThan = sum(avail,2)>=15 & sum(avail(:,end-2:end),2)>0;
finalSelect= SU2_acrossSes  & availMoreThan;


rng(1234)
switch nVS
    case 1
groundTruthStimPattern_linear=reshape(repmat([1:16],50,1)',[],1);
nSuper=16
totaltrails=800
    case 5
groundTruthStimPattern_linear=reshape(repmat([1:100],60,1)',[],1);
nSuper=100
totaltrails=6000
case 3
groundTruthStimPattern_linear=reshape(repmat([1:81],30,1)',[],1);
nSuper=81
totaltrails=2430
    case 2
groundTruthStimPattern_linear=reshape(repmat([1:30],90,1)',[],1);
nSuper=30      
totaltrails=2700
end
train_trials=arrayfun(@(x) sri(totaltrails,x,1:floor(totaltrails*6/7)),1:7,'UniformOutput',false);
train_trials=cat(2,train_trials{:});

valid_trails=setdiff(1:(totaltrails*7),train_trials);

ErrorArrayAnimal=cell(1,5);
decoderInfoAnimal=cell(3,5,15);

for ani = 1:5

    selUid=find(finalSelect & whichAnimalPerUnit==ani & umapWorstDistance<prctile(umapWorstDistance(whichAnimalPerUnit==ani),umapThres));
    decoderInfoAnimal{1,ani,1}=selUid;
    fc=VS_Store(selUid,:);

for day = 8:15 % 
    daysSelected=day-7:day-1;

dayID = repmat(daysSelected,totaltrails,1);
dayID=dayID(:);




    X=[];
    Xnew=[];
    for dayS=1:numel(daysSelected)
        X=[X cell2mat(cellfun(@(x) reshape(x,[],1)',fc(:,daysSelected(dayS)),'UniformOutput', false))];
    end
    YY = repmat(groundTruthStimPattern_linear,numel(daysSelected),1);
    XXX=X(:,valid_trails);YYY=YY(valid_trails);X=X(:,train_trials);YY=YY(train_trials);


       [Mdl] = fitcdiscr(fillmissingHanlin(X)',YY,...
    'OptimizeHyperparameters',{'Gamma'},...
    'HyperparameterOptimizationOptions',struct('ShowPlots',false, ...
    'UseParallel',false,'KFold',5,...
    'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30,'Repartition',true));   
       Mdl.Gamma=Mdl.HyperparameterOptimizationResults.XAtMinObjective.Gamma;

decoderInfoAnimal{2,ani,day}=Mdl.Gamma;

    Xnew=XXX;
    badrows=sum(isnan(Xnew),2)>0;
    Xnew=Xnew(~badrows,:);

    Xnew=fillmissingHanlin([Xnew X(~badrows,:)]);
    Xnew=Xnew(:,1:numel(YYY));
    maxI=predict(Mdl,Xnew')';
    decodedAngle = (360*(maxI-1)/nSuper)';

    ydegreeNew=(YYY-1)*22.5;
    error = abs(decodedAngle - ydegreeNew);
    error(error>180)=360-error(error>180);
    ErrorArray{4,day}=error;
    ErrorArray{5,day}=maxI'==YYY;

    [day mean(error) mean(maxI'==YYY)]

    
    Xnew=cell2mat(cellfun(@(x) reshape(x,[],1)', fc(:,day),'UniformOutput', false)); 
    badrows=sum(isnan(Xnew),2)>0;
    Xnew=Xnew(~badrows,:);
    Mdl = fitcdiscr(fillmissingHanlin(X(~badrows,:))',YY,'Gamma',Mdl.Gamma);
    Xnew=fillmissingHanlin([Xnew X(~badrows,:)]);
    Xnew=Xnew(:,1:numel(groundTruthStimPattern_linear));
    maxI=predict(Mdl,Xnew')';
    decodedAngle = (360*(maxI-1)/nSuper)';
    ydegreeNew= (groundTruthStimPattern_linear-1)*22.5; % ground truth new data
    error = abs(decodedAngle - ydegreeNew);
    error(error>180)=360-error(error>180);
    ErrorArray{1,day}=error;
    ErrorArray{2,day}=maxI'==groundTruthStimPattern_linear;

    [day mean(error) mean(maxI'==groundTruthStimPattern_linear)]

end
ErrorArrayAnimal{ani}=ErrorArray;
end
save(['ErrorArrayAnimal_ratesOnly_Slide.mat'],'ErrorArrayAnimal','decoderInfoAnimal')
%% ploting
try
    LMEtable;
catch
    LMEtable=table;
end
foldOverChance=1; 
load('nonNoise_whichAnimalLayerShankID7.mat')


switch nVS
    case 1
nSuper=16;
    case 5
nSuper=100;
case 3
nSuper=81;
    case 2
nSuper=30;      
end
if foldOverChance==1
chance=1/nSuper;
else
    chance=1;
end


load('ErrorArrayAnimal_ratesOnly_Slide.mat')
driftRateCP=cell(5,8);
for ani =1:5
driftRateCP(ani,:)=cellfun(@(x) mean(reshape(x(end-floor(numel(x)/nSuper)*nSuper+1:end),nSuper,[]),2), ErrorArrayAnimal{ani}(3-foldOverChance,8:15),UniformOutput=false);
end
driftRateCP=cell2mat(driftRateCP);

driftRateWithin=cell(5,8);
for ani =1:5
driftRateWithin(ani,:)=cellfun(@(x) mean(reshape(x(end-floor(numel(x)/nSuper)*nSuper+1:end),nSuper,[]),2), ErrorArrayAnimal{ani}(6-foldOverChance,8:15),UniformOutput=false);
end

driftRateWithin=cell2mat(driftRateWithin);

load('ErrorArrayAnimal_ratesOnly.mat')
    driftRate=cell(5,8);
    for ani =1:5
    driftRate(ani,:)=cellfun(@(x) mean(reshape(x(end-floor(numel(x)/nSuper)*nSuper+1:end),nSuper,[]),2), ErrorArrayAnimal{ani}(3-foldOverChance,8:15),UniformOutput=false);
    end
   
    driftRate=cell2mat(driftRate);
 

load('ErrorArrayAnimal_ratesOnly_Shuf.mat')
 driftRateShuf=cell(5,8);
    for ani =1:5
    driftRateShuf(ani,:)=cellfun(@(x) mean(reshape(x(end-floor(numel(x)/nSuper)*nSuper+1:end),nSuper,[]),2), ErrorArrayAnimal{ani}(3-foldOverChance,8:15),UniformOutput=false);
    end
   
    driftRateShuf=cell2mat(driftRateShuf);



    driftRate=100*(1-driftRate);
    driftRateShuf=100*(1-driftRateShuf);
    driftRateCP=100*(1-driftRateCP);
    driftRateWithin=100*(1-driftRateWithin);

h=figure();
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 7.5;
h.Position(4)         = 11;

chance=(1-chance)*100;
plot([1 8],[chance chance],'k--')
hold on
    hanlinErrorBar(driftRateShuf,'m')

    hanlinErrorBar(driftRate,'b')
     hold on;
    hanlinErrorBar(driftRateWithin,'#FFA500')
    hold on
    hanlinErrorBar(driftRateCP,'r')
    hold on;
    xlabel('interval/day')
    
    ylabel('decoding error %'); set(gca,'YDir','normal');
    ylim([0 100])
 yl=ylim;

yl=ylim;
xlim([0.5 8.5])
xticks([1:2:9])
xticklabels(num2cellstr(1:2:9))                    
    legend('chance','shuffle','fixed','within day','sliding','EdgeColor','none','Color','none','Location','northwest')
    set(gca,'FontSize',16)
    title(['decode ' stimNames{nVS}],'FontWeight','normal')
if nVS==5
  title(['decode ' '100 images'],'FontWeight','normal','FontSize',12)
end
set(findall(h,'-property','FontSize'),'FontSize',12)
set(findall(h,'Type','Legend','-property','FontSize'),'FontSize',11)

box off
lme=hanlinLME([{driftRateCP} {driftRate} ],1);

str1=sprintf('p(\\color{red}same \\color{blue}slopes\\color{black})=%.2g',double(lme{4}.Coefficients(4,6)) );
str2=sprintf('p(\\color{red}slope \\color{black}= 0)=%.2g',double(lme{4}.Coefficients(2,6)) );
text(0.75,yl(2)*0.25,str2,'FontSize',11)
text(0.75,yl(2)*0.18,str1,'FontSize',11)

entry=dataset2table(lme{4}.Coefficients(4,:));
entry.N=size(driftRate,1);entry.Unit={'mice-pattern pairs'};entry.Fig={['d-' stimNames{nVS} '-fssc']};
LMEtable(end+1,:)=entry;

entry=dataset2table(lme{4}.Coefficients(2,:));
entry.N=size(driftRate,1);entry.Unit={'mice-pattern pairs'};entry.Fig={['d-' stimNames{nVS} '-fssc']};
LMEtable(end+1,:)=entry;

lme=hanlinLME([{driftRate} {driftRateShuf} ],1);

str2=sprintf('p(\\color{blue}same \\color{magenta}slopes\\color{black})=%.2g',double(lme{4}.Coefficients(4,6)));
if double(lme{4}.Coefficients(3,6))>realmin
    str1=sprintf('p(\\color{blue}same \\color{magenta}intercept\\color{black})=%.2g',double(lme{4}.Coefficients(3,6)) );
else
str1=sprintf('p(\\color{blue}same \\color{magenta}intercept\\color{black})<%.2g',realmin);
end
text(0.75,yl(2)*0.11,str2,'FontSize',11)
text(0.75,yl(2)*0.04,str1,'FontSize',11)

entry=dataset2table(lme{4}.Coefficients(3,:));
entry.N=size(driftRate,1);entry.Unit={'mice-pattern pairs'};entry.Fig={['d-' stimNames{nVS} '-fssc']}
LMEtable(end+1,:)=entry;

entry=dataset2table(lme{4}.Coefficients(4,:));
entry.N=size(driftRate,1);entry.Unit={'mice-pattern pairs'};entry.Fig={['d-' stimNames{nVS} '-fssc']}
LMEtable(end+1,:)=entry;

if nVS==1
ylim([-10 100])
end
if nVS==2
ylim([10 100])
end
if nVS==3
ylim([60 100])
end
storeFolder = pwd;
if nVS==5
saveas(h,[storeFolder filesep entry.Fig{1} '.fig'])
exportgraphics(h,[storeFolder filesep entry.Fig{1} '.pdf'],'ContentType','vector')
else
close(h);
end
%%
h=figure('Position',[819   408   368 377]);

    hanlinErrorBar(driftRate,'b')
     hold on;
    hanlinErrorBar(driftRateWithin,'#FFA500')
    hold on
    hanlinErrorBar(driftRateCP,'r')
    hold on;
    xlabel('interval/day')
    ylabel('decoding error %'); set(gca,'YDir','normal');
    ylim([0 100])
    if nVS==3
ylim([50 100])
end
 yl=ylim;
yl=ylim;
xlim([0.5 8.5])
xticks([1:2:9])
xticklabels(num2cellstr(1:2:9))                    
    legend('across day fixed','within day','across day slide','EdgeColor','none','Color','none','Location','northwest')
    set(gca,'FontSize',16)
    title(['decode ' stimNames{nVS}],'FontWeight','normal')
if nVS==5
  title(['decode ' '100 images'],'FontWeight','normal','FontSize',12)
end


box off
lme=hanlinLME([{driftRateCP} {driftRate} ],1);

str1=sprintf('p(\\color{red}same \\color{blue}slopes\\color{black})=%.2g',double(lme{4}.Coefficients(4,6)) );
str2=sprintf('p(\\color{red}slope \\color{black}= 0)=%.2g',double(lme{4}.Coefficients(2,6)) );


text(0.75,yl(1)+range(yl)*0.15,str2,'FontSize',15)
text(0.75,yl(1)+range(yl)*0.05,str1,'FontSize',15)


  titlestr=['LDA_decode_Rate_VS_Slide' stimNames{nVS}];
saveas(h,titlestr,'fig');
saveas(h,titlestr,'png');




if foldOverChance==1
chance=1/nSuper;
else
    chance=1;
end
% plot shuffled ------------------


h=figure('Position',[819   408   368 377]);
chance=(1-chance)*100;
plot([1 8],[chance chance],'k--')
hold on
    hanlinErrorBar(driftRateShuf,'m')
     hold on;
    hanlinErrorBar(driftRate,'b')
  
    xlabel('interval/day')
    ylabel('decoding error %'); set(gca,'YDir','normal');
    ylim([0 100])
    if nVS==3
ylim([50 100])
end
yl=ylim;

yl=ylim;
xlim([0.5 8.5])
xticks([1:2:9])
xticklabels(num2cellstr(1:2:9))                     
    legend('chance','shuffle','fixed','EdgeColor','none','Color','none','Location','northwest')
    set(gca,'FontSize',16)
lme=hanlinLME([{driftRate} {driftRateShuf} ],1);

str1=sprintf('p(\\color{blue}same \\color{magenta}intercept\\color{black})=%.2g',double(lme{4}.Coefficients(3,6)) );
str2=sprintf('p(\\color{blue}same \\color{magenta}slope\\color{black})=%.2g',double(lme{4}.Coefficients(4,6)) );


box off
if nVS==3
ylim([50 100])
end
text(0.75,yl(1)+range(yl)*0.15,str2,'FontSize',15)
text(0.75,yl(1)+range(yl)*0.05,str1,'FontSize',15)

        titlestr=['LDA_decode_Rate_VS_Shuf' stimNames{nVS}];
saveas(h,titlestr,'fig');
saveas(h,titlestr,'png');
end
