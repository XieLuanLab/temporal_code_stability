%% OOO rates only LDA based decoding analysis
basePath = pwd;
addpath(genpath(basePath));
for nVS=[1 2 3 5]
    nthComp=1:3;
    NumComp=3;
    daysSelected=1:7;

    usePCA=false;
    perDayCVcomp=true;
    ZeroMean=0;
    DivideSum=1;
    kernalWidth=24;
    kfold=5;

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

    availMoreThan = sum(avail,2)>=12 & all(avail(:,(max(daysSelected)+1):end),2);

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


    nTrialsPerDay=totaltrails/nSuper;
    sample1outofNtrialsPerDay=numel(daysSelected);
    nTrialsPerDayHoldOut = floor(nTrialsPerDay/sample1outofNtrialsPerDay);
    testTrialIndex=[];
    for  day=1:numel(daysSelected)
        for batch =1:nTrialsPerDayHoldOut
            testTrialIndex(end+1)=(day-1)*nTrialsPerDay+randi(sample1outofNtrialsPerDay)+(batch-1)*sample1outofNtrialsPerDay;
        end
    end

    valid_trails=arrayfun(@(x) sri(nSuper,x,1:nSuper),testTrialIndex,'UniformOutput',false);
    valid_trails=cat(2,valid_trails{:});
    train_trials=setdiff(1:(totaltrails*numel(daysSelected)),valid_trails);



    groundTruthStimPattern_linear_test=groundTruthStimPattern_linear;
    gTSPLA=1:numel(groundTruthStimPattern_linear_test);
    if nVS==3
        valid=[21:2:27 39:2:45 57:2:63 75:2:81];
        validTrials=bsxfun(@plus,repmat([0:81:(numel(daysSelected)*2430-81)]',1,16),valid)';
        [LA,LB]=ismember(train_trials,validTrials(:));
        train_trials=train_trials(LA);
        [LA,LB]=ismember(valid_trails,validTrials(:));
        valid_trails=valid_trails(LA);
        [gTSPLA,LB]=ismember(groundTruthStimPattern_linear,valid);
        groundTruthStimPattern_linear_test=groundTruthStimPattern_linear(gTSPLA);
    end


    cvOptions=cvpartition("CustomPartition",dayID(train_trials)-min(dayID)+1);


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


    costMatrix=[];
    switch nVS
        case 1
            costMatrix=pdist2([0:22.5:337.5]',[0:22.5:337.5]');
            costMatrix(costMatrix>180)=360-costMatrix(costMatrix>180);
        case 2
            costMatrix=[];
        case 3
            costMatrix=AngleErrorMatrix;
            costMatrix=costMatrix(valid,valid);
        case 5
            costMatrix=[];
    end

    %% firing rate based decoding
    for ani=[1:5]
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
        [Mdl] = fitcdiscr(data_X,YY,'Cost',costMatrix,...
            'OptimizeHyperparameters',{'Gamma'},...
            'HyperparameterOptimizationOptions',struct('ShowPlots',false, ...
            'UseParallel',false,'CVPartition',cvOptions,...
            'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30));
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
            Mdl = fitcdiscr(data_X,YY,'Cost',costMatrix,'Gamma',Mdl.Gamma);
            [maxI,IdxProb]=predict(Mdl,Xnew');maxI=maxI(:);maxI=maxI';
            if day~=7 && nVS==3
                maxI=maxI(gTSPLA);
                IdxProb= IdxProb(gTSPLA,:);
            end

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
    save(['EAA_ratesOnly.mat'],'ErrorArrayAnimal','decoderInfoAnimal')
    %% with temporal components decoding
    score7=[];
    for ani = 1:5
        ani
        filename=['LDA_' num2str(kernalWidth*2) 'ms_' num2str(kfold) 'pD' num2str(perDayCVcomp) 'CV_a_' num2str(ani) 'nVS_' num2str(nVS) '_ZeroMean_' num2str(ZeroMean) '_DivideSum_' num2str(DivideSum) '_PCA' num2str(usePCA) '.mat'];
        load(filename,'LDA_7score_Store')
        LDA_7score_Store=fillEmptyAsNaN(LDA_7score_Store);
        LDA_7score_Store=cellfun(@(x) x(:,:,nthComp),LDA_7score_Store,'UniformOutput',false);
        score7=[score7;LDA_7score_Store];
    end
    LDA_7score_Store=[];
    score7(avail==0)={nan(size(score7{1}))};

    decoderInfoAnimal=cell(3,5);
    ErrorArrayAnimal=cell(1,5);
    for ani=[1:5]

        selUid=find(finalSelect & whichAnimalPerUnit==ani & umapWorstDistance<prctile(umapWorstDistance(whichAnimalPerUnit==ani),umapThres));

        ft=score7(selUid,:);
        ft=cellfun(@real,ft,'UniformOutput',false);
        fc=VS_Store(selUid,:);


        XX=[];
        Xnew=[];
        for comp=1:NumComp
            X=[];
            for dayS=1:numel(daysSelected)
                X=[X cell2mat(cellfun(@(x) reshape(x(:,:,comp),[],1)',ft(:,daysSelected(dayS)),'UniformOutput', false))];
            end
            XX=[XX;X];
        end
        selCid=repmat(selUid(:),NumComp,1);

        X=[];
        for dayS=1:numel(daysSelected)
            X=[X cell2mat(cellfun(@(x) reshape(x,[],1)',fc(:,daysSelected(dayS)),'UniformOutput', false))];
        end
        XX=[XX;fillmissingHanlin(X,true)];
        selCid=[selCid;selUid(:)];

        filledXX=fillmissingHanlin(XX,false);

        badrowsXX=var(filledXX,[],2)<0.001 | isnan(var(filledXX,[],2));
        XX=XX(~badrowsXX,:);
        selCid=selCid(~badrowsXX);
        decoderInfoAnimal{1,ani}=selCid;

        YY = repmat(groundTruthStimPattern_linear,numel(daysSelected),1);

        XXX=XX(:,valid_trails); YYY=YY(valid_trails); XX=XX(:,train_trials);YY=YY(train_trials);


        data_X=fillmissingHanlin(XX,false)';


        [Mdl] = fitcdiscr(data_X,YY,'Cost',costMatrix,...
            'OptimizeHyperparameters',{'Gamma'},...
            'HyperparameterOptimizationOptions',struct('ShowPlots',false, ...
            'UseParallel',false,'CVPartition',cvOptions,...
            'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30));
        Mdl.Gamma=Mdl.HyperparameterOptimizationResults.XAtMinObjective.Gamma;

        decoderInfoAnimal{2,ani}=Mdl.Gamma;


        for day =7:15
            day
            XXnew=[];
            for comp=1:NumComp
                Xnew=cell2mat(cellfun(@(x) reshape(x(:,:,comp),[],1)',ft(:,day),'UniformOutput', false));
                XXnew=[XXnew;Xnew];
            end
            XXnew=[XXnew; cell2mat(cellfun(@(x) reshape(x,[],1)',fc(:,day),'UniformOutput', false))];
            XXnew=XXnew(~badrowsXX,:);

            if day==7
                XXnew=XXX;
            end

            badrows=sum(isnan(XXnew),2)<0;% all good rows.
            XXnew=XXnew(~badrows,:);
            data_X=fillmissingHanlin(XX(~badrows,:),false)';%this shouldn't do anything
            XXnew=fillmissingHanlin([XXnew XX(~badrows,:)],false);%this shouldn't do anything

            if day==7
                XXnew=XXnew(:,1:numel(YYY));
            else
                XXnew=XXnew(:,1:numel(groundTruthStimPattern_linear));
            end

            Mdl = fitcdiscr(data_X,YY,'Gamma',Mdl.Gamma,'Cost',costMatrix);

            [maxI,IdxProb]=predict(Mdl,XXnew');maxI=maxI(:);maxI=maxI';

            if day~=7 && nVS==3
                maxI=maxI(gTSPLA);
                IdxProb= IdxProb(gTSPLA,:);
            end

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


            if nVS==3 && day==7 %
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
    save(['EAA_CompRates.mat'],'ErrorArrayAnimal','decoderInfoAnimal')
%% decode plot 
CPstr='EAA_CompRates'
FRstr='EAA_ratesOnly'
global xaxisshaded
xaxisshaded=2:9;
load('nonNoise_whichAnimalLayerShankID7.mat')
 stimNames = {'DG','SG','RFG','','IJO'};
cd([basePath filesep stimNames{nVS}])

switch nVS
    case 1
nSuper=16;
foldOverChance=2; % 
    case 5
nSuper=100;
foldOverChance=1; % 
case 3
nSuper=25;
foldOverChance=2; % 
    case 2
nSuper=30;      
foldOverChance=1; % 
end

if foldOverChance==1
chance=1/nSuper;
else
    chance=1;
end

load([ pwd filesep CPstr '.mat'])

driftRateCP=cell(5,9)
for ani =1:5  
driftRateCP(ani,:)=cellfun(@(x) mean(reshape(x(end-floor(numel(x)/nSuper)*nSuper+1:end),nSuper,[]),2), ErrorArrayAnimal{ani}(3-foldOverChance,7:15),UniformOutput=false);
end

driftRateCP=double(cell2mat(driftRateCP)/chance);

load([ pwd filesep FRstr '.mat'])

    driftRate=cell(5,9)
    for ani =1:5
    driftRate(ani,:)=cellfun(@(x) mean(reshape(x(end-floor(numel(x)/nSuper)*nSuper+1:end),nSuper,[]),2), ErrorArrayAnimal{ani}(3-foldOverChance,7:15),UniformOutput=false);
    end

    driftRate=double(cell2mat(driftRate)/chance);

h=figure();
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 4.0*2*1.05;
h.Position(4)         = 4.0*2;
h.PaperPositionMode   = 'auto';

hanlinErrorBar(driftRateCP(:,min(xaxisshaded):end),'r')
hold on;
hanlinErrorBar(driftRate(:,min(xaxisshaded):end),'b')
hold on;

H=shadedErrorBar([0.75 1.25],repmat(nanmean(driftRateCP(:,1)),1,2),sem(driftRateCP(:,1)),'lineProps',{'lineStyle','none','Color','r'});
H.edge(2).Visible='off';
H.edge(1).Visible='off';
H=shadedErrorBar([0.75 1.25],repmat(nanmean(driftRate(:,1)),1,2  ),    sem(driftRate(:,1)),'lineProps',{'lineStyle','none','Color','b'});
H.edge(2).Visible='off';
H.edge(1).Visible='off';
scatter(1,nanmean(driftRateCP(:,1)),25,'r','Marker','o','MarkerFaceColor','white','MarkerEdgeColor','red','LineWidth',2)
scatter(1,nanmean(driftRate(:,1)),25,'b','Marker','o','MarkerFaceColor','white','MarkerEdgeColor','blue','LineWidth',2)

    xlabel('interval (days)')
    ylabel('accuracy (fold over chance)')
    if (nVS==1 || nVS==3) && foldOverChance==2
    ylabel('error (degree)'); set(gca,'YDir','reverse');ylim([7 50])
    end

if nVS==5
     ylim([40 85])
end
if nVS==2
    ylim([11 22])
end
if nVS==3 && foldOverChance==2
    ylim([7 19])
end
yl=ylim;
xlim([0.5 9.5])
xticks([1:2:9])
xticklabels(num2cellstr(0:2:8))                     
    leg=legend('components+rate','rate','Color','none','EdgeColor','none');
    leg.Position=[0.3529 0.7865 0.5570 0.1426];
    set(gca,'FontSize',11)
    title(['decode ' stimNames{nVS}],'FontWeight','normal','FontSize',12)
if nVS==5
  title(['decode ' 'NI'],'FontWeight','normal','FontSize',12)
end
set(findall(h,'-property','FontSize'),'FontSize',12)
set(findall(h,'Type','Legend','-property','FontSize'),'FontSize',11)
    lme=hanlinLME([{driftRate(:,[min(xaxisshaded):end])} {driftRateCP(:,[min(xaxisshaded):end])}],1);
    if nVS==1 && foldOverChance==2
ylim([7 50])
    end
str1=sprintf('p(same slopes)=%.2g',double(lme{4}.Coefficients(4,6)));
str2=sprintf('p(same intercepts)=%.2g',double(lme{4}.Coefficients(3,6)));

if nVS==2 || (nVS==1 && foldOverChance~=2) || nVS==5
text(2,range(yl)*0.14+yl(1),str2,'FontSize',11,'Color','k')
text(2,range(yl)*0.055+yl(1),str1,'FontSize',11,'Color','k')
elseif (nVS==3 && foldOverChance==2)
text(2,17.5,str2,'FontSize',11,'Color','k')
text(2,18.5,str1,'FontSize',11,'Color','k')
else
text(2,43,str2,'FontSize',11,'Color','k')
text(2,47,str1,'FontSize',11,'Color','k')
end
xaxisshaded=[]
%% 
entry.Fig={[stimNames{nVS} '-deocde-FR-CP']};
storeFolder = pwd;
saveas(h,[storeFolder filesep entry.Fig{1} '.fig'])
exportgraphics(h,[storeFolder filesep entry.Fig{1} '.pdf'],'ContentType','vector')
%% 
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
% put in vector form
a = [[cosd(AzEl1(1)) sind(AzEl1(1))] * cosd(AzEl1(2)) sind(AzEl1(2))];
b = [[cosd(AzEl2(1)) sind(AzEl2(1))] * cosd(AzEl2(2)) sind(AzEl2(2))];
AngleErrorMatrix(i,j) = acosd(dot(a,b));
        end
        end
end

if nVS==1 && foldOverChance==2
    driftRate=180-driftRate;
    driftRateCP=180-driftRateCP;     
end
if nVS==3 && foldOverChance==2
try
    valid=[21:2:27 39:2:45 57:2:63 75:2:81]
    driftRate=max(AngleErrorMatrix(valid,valid),[],'all')-driftRate;
    driftRateCP=max(AngleErrorMatrix(valid,valid),[],'all')-driftRateCP; 
catch
    driftRate=max(AngleErrorMatrix(:))-driftRate;
    driftRateCP=max(AngleErrorMatrix(:))-driftRateCP;     
end
end
driftRateMG=(driftRate(:,1)-driftRate(:,1:end))./(driftRate(:,1)+driftRate(:,1:end));
driftRateCPMG=(driftRateCP(:,1)-driftRateCP(:,1:end))./(driftRateCP(:,1)+driftRateCP(:,1:end));
%% 
h=figure();
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 4.0*2*1.05;
h.Position(4)         = 4.0*2;
h.PaperPositionMode   = 'auto';

hanlinErrorBar(driftRateCPMG(:,[1 2:9]),'r')
hold on;
hanlinErrorBar(driftRateMG(:,[1 2:9]),'b')

H=shadedErrorBar([0.75 1.25],repmat(nanmean(driftRateCPMG(:,1)),1,2),sem(driftRateCPMG(:,1)),'lineProps',{'lineStyle','none','Color','r'});
H.edge(2).Visible='off';
H.edge(1).Visible='off';
H=shadedErrorBar([0.75 1.25],repmat(nanmean(driftRateMG(:,1)),1,2  ),    sem(driftRateMG(:,1)),'lineProps',{'lineStyle','none','Color','b'});
H.edge(2).Visible='off';
H.edge(1).Visible='off';
scatter(1,nanmean(driftRateCPMG(:,1)),25,'r','Marker','o','MarkerFaceColor','white','MarkerEdgeColor','red','LineWidth',2)
scatter(1,nanmean(driftRateMG(:,1)),25,'b','Marker','o','MarkerFaceColor','white','MarkerEdgeColor','blue','LineWidth',2)

hold on;
xlabel('interval (days)')
ylabel('decoding drift index')
switch nVS
    case 1
ylim([0 0.1])
 case 2
     ylim([0 0.18])
     case 3
     ylim([0 0.065])
     case 5
     ylim([0 0.21])
end
yl=ylim;

xlim([0.5 9.5])
xticks([1:9])
xticklabels(num2cellstr(0:8))                     
leg=legend('components+rate','rate','Color','none','EdgeColor','none','Location','northwest');
leg.Position=[0.2093    0.7948    0.5830    0.1390];
    set(gca,'FontSize',12)
    title(['decode ' stimNames{nVS}],'FontWeight','normal')
if nVS==5
  title(['decode ' 'NI'],'FontWeight','normal')
end
    lme=hanlinLME([{driftRateMG(:,[1 2:9])} {driftRateCPMG(:,[1 2:9])}],0);

    set(findall(h,'-property','FontSize'),'FontSize',12)
set(findall(h,'Type','Legend','-property','FontSize'),'FontSize',11)

str1=sprintf('p(same slopes)=%.2g',string(double(lme{4}.Coefficients(4,6))));
str2=sprintf('p(same intercepts)=%.2g',string(double(lme{4}.Coefficients(3,6))));
        text(2.5,yl(2)*0.05,str1,'FontSize',11,'Color','k')

entry.Fig={[stimNames{nVS} '-deocdeRI-FR-CP']};
storeFolder = pwd;
saveas(h,[storeFolder filesep entry.Fig{1} '.fig'])
exportgraphics(h,[storeFolder filesep entry.Fig{1} '.pdf'],'ContentType','vector')

end
