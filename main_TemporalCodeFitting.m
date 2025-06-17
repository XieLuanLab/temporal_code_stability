mkdir('results')
% main code start
% VS_Store = cellfun(@(x) x(:,:,1:2:end)+x(:,:,2:2:end),VS_Store ,'Uni',0); % downsample 1ms 

load('VS_Store.mat')
nU=size(VS_Store,1);
nD=size(VS_Store,2);
LDA_15Axis_Store=cell(nU,1);
LDA_15CV_Store=nan(nU,1);
LDA_15gamma_Store=nan(nU,1);
LDA_15lambda_Store=cell(nU,1);
LDA_15score_Store=cell(nU,nD);
LDA_decodePower_trialS=nan(nU,1);

trialN=50;% 50 trials
stimN=16;% 16 grating directions

u = 1; % select this example unit
avail = true(1,15); % this unit showed up for all 15 days.
rng(1234)
label=reshape(repmat(1:stimN,trialN*sum(avail(u,:)),1)',[],1);
temp=cellfun(@(x) permute(x,[2 3 1]),VS_Store(u,:),'UniformOutput',false);
tempSum = cellfun(@(x) squeeze(sum(x,2)),temp,'UniformOutput',false); % total firing count
temp=cellfun(@(x) smoothdata(x,2,'movmean',[24 24]),temp,'UniformOutput',false); % boxcar filter 49 ms
temp=cellfun(@(x) x(:, floor(1:2:end),:),temp,'UniformOutput',false); % downsample by half
firingTimeCourse=cat(3,temp{avail(u,:)}); % remove day axis , as if one day
firingTimeCourse=permute(firingTimeCourse,[2,1,3]);
X0 = firingTimeCourse(:,:); % X0 of shape timeBin X trials

Xs=sum(X0,1);
Xs(Xs==0)=1;
X0=bsxfun(@rdivide,X0,Xs);

dayID = repmat(1:(numel(label)/trialN/stimN),trialN*stimN,1);
dayID=dayID(:);
cvOptions=cvpartition("CustomPartition",dayID-min(dayID)+1);

params = hyperparameters('fitcdiscr',X0',label);
params(1)=[];

badRows=any(cell2mat(arrayfun(@(par) max(X0(:,training(cvOptions,par)),[],2)==0,1:cvOptions.NumTestSets,'UniformOutput',false)),2); 

try
   [Mdl] = fitcdiscr(X0(~badRows,:)',label,...
    'OptimizeHyperparameters',params,...
    'HyperparameterOptimizationOptions',struct('ShowPlots',false, ...
    'UseParallel',false,'CVPartition',cvOptions,...
    'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30,'Repartition',false));
    [Mdl ] = fitcdiscr(X0(~badRows,:)',label,'DiscrimType','linear','Gamma',Mdl.HyperparameterOptimizationResults.XAtMinObjective.Gamma);
catch % fallback to diagonal linear if failed
    [Mdl ] = fitcdiscr(X0',label,'DiscrimType','diaglinear');
end

CVMdl=[];
CVMdl.kfoldLoss=lda_kfold_cv(X0, label, Mdl.Gamma,5,7);
LDA_15gamma_Store(u)=Mdl.Gamma;

try
    [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma);
catch % fallback to diagonal linear if failed
    [W, LAMBDA] = eig(Mdl.BetweenSigma, diag(Mdl.Sigma));
end

lambda = diag(LAMBDA);
[lambda, SortOrder] = sort(lambda, 'descend','MissingPlacement','last'); % sort LDA components
W = W(:, SortOrder);W=W(:,1:min((stimN-1),size(W,2)));if any(~badRows);W=interp1(find(~badRows),W,[1:size(X0,1)]','linear',0);end
 % get LDA weight for timeBins (components)
Y = (X0')*W;  Z=reshape(Y,stimN,trialN,[],(stimN-1));% get LDA scores after projecting along components N-Class - 1

LDA_15Axis_Store{u}=W;
LDA_15CV_Store(u)=1-CVMdl.kfoldLoss;

dlist=find(avail(u,:));
for u1=1:numel(dlist) % get per day LDA scores
    LDA_15score_Store{u,dlist(u1)}=squeeze(Z(:,:,u1,:));
    LDA_15score_Store{u,dlist(u1)}=LDA_15score_Store{u,dlist(u1)}(:,:,1:7); % only save 7 components
end
firingRates=cat(3,tempSum{dlist}); % get firing rate based decoding error
X=firingRates(:);
try
    [Mdl ] = fitcdiscr(X,label,'DiscrimType','linear');
catch
    [Mdl ] = fitcdiscr(X,label,'DiscrimType','diaglinear');
end
CVMdl=[];
CVMdl.kfoldLoss=lda_kfold_cv(X', label, 0, 5,7);
LDA_decodePower_trialS(u)=1-CVMdl.kfoldLoss;

%% plot Fig4a
VS_Store_Avg = cellfun(@(x) squeeze(mean(x(:,:,1:500),1)),VS_Store ,'Uni',0);
h=figure();
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 26;
h.Position(4)         = 9.6;
h.PaperPositionMode   = 'auto';
tl=tiledlayout(4,4,'TileSpacing','tight','Padding','compact')

colorList=gray(15);
clist=hsv(16);
clist=[clist(1:2:end,:);clist(2:2:16,:)];
temp=cellfun(@(x) x,VS_Store_Avg(u,:),'UniformOutput',false);
temp=cellfun(@(x) smoothdata(x,2,'movmean',[24 24]),temp,'UniformOutput',false);
ylMax=max(cellfun(@(x) max(x(:)),temp))*1000;
for stim=1:16
    nexttile
    for d=1:15
        plot(temp{d}(stim,:)*1000,'-','Color',[colorList(d,:) 0.7]);
        hold on
    end
    ylim([0 ylMax])
    xlim([0 500])
    yticks([])
    xticks([])
    text(500,ylMax,[num2str((stim-1)*22.5) '°'],'FontSize',12,'VerticalAlignment','top','HorizontalAlignment','right')
    if stim==1
        yticks('auto')
    end
    if stim==16
        xticks('auto')
    end
    box off
end



h.Position(3)         = h.Position(3)/2;
h.Position(4)         = h.Position(4)/2;
c=colorbar;
colormap(colorList);
c.Ticks=[1/15:1/15:1]-1/30
c.Position=[0.8332 0.1195 0.02 0.84]
c.TickLabels=num2cellstr(1:15)
ylabel(c,'day');
tl.Position=[0.0650    0.1195    0.7550    0.8430]
set(findall(h,'-property','FontSize'),'FontSize',6)

exportgraphics(h,'results/Fig4a-DG-PSTH.pdf','ContentType','vector')
close(h)


%% plot Fig4b
h=figure();
 h.Units               = 'centimeters';
 h.Position(1)         = 3;
 h.Position(2)         = 3;
 h.Position(3)         = 5.6*2;
 h.Position(4)         = 4.2*2;
 h.PaperPositionMode   = 'auto';

 plot([1 500],[0 0],'k:','lineWidth',1)
 hold on
 newcolors = {'#F00','#00F','#A0F'};
for i =1:3
    thisCP=LDA_15Axis_Store{u}(1:250,i);
    if max(thisCP)<abs(min(thisCP))
    thisCP=-thisCP;
    end
plot(1:2:500,thisCP,'lineWidth',1.5,'Color',newcolors{i})    
hold on
end
ylim([-30 30])
set(gca,'FontSize',12)
 xlabel('time ms')
 ylabel('weight')
l=legend('','CP1','CP2','CP3','EdgeColor','none','Location','north');
l.Position=[0.4332    0.7130    0.1951    0.2183];
xlim([0 520])


set(findall(h,'-property','FontSize'),'FontSize',14)
box off
exportgraphics(h,'results/Fig4b-DG-LDA-CP.pdf','ContentType','vector')
close(h);

%% plot Fig4c/left
figure
clist=hsv(16);
clist=[clist(1:2:end,:);clist(2:2:16,:)];
clist(3,:) = [121 130 56]/255;
clist(4,:) = [1 121 111]/255;

for i = [16:-1:1]

    temp=cellfun(@(x) sum(x(:,i,:),[2 3]),VS_Store(u,:),'UniformOutput',0);
    temp=cat(2,temp{:})'/0.52; % firing rate of size 15d x 50 trial

    ascomb=sqrt(temp+3/8)*2;
    ascomb=mean(ascomb,2); % average across trials

    pd=fitdist(ascomb,'Normal');
    xpdf=linspace(pd.mu-pd.sigma,pd.mu+pd.sigma,30);
    hold on;area(xpdf,pdf(pd,xpdf)/max(pdf(pd,xpdf)),'FaceColor',clist(i,:),'FaceAlpha',1/8,'EdgeColor',clist(i,:),'EdgeAlpha',0.8,'LineWidth',2)
end
ylabel('Normalized Probability')
xlabel('Transformed Firing Rate (Hz)')
set(gca,'FontSize',18)
ylim([0.606 1])
yticks([0.7 0.8 0.9 1])
xlim([3.6 6.3])
%xlim([2.7 4.6])
exportgraphics(gcf,'results/Fig4c-left.pdf','ContentType','vector')

saveas(gcf,'results/Fig4c-left.png')
close(gcf)

%% plot Fig4c/right
figure('Position',[679     1   644   420])
for i = [16:-1:1]
    temp=cellfun(@(x) x(i,:,1:3),LDA_15score_Store(u,:),'UniformOutput',0);
    temp=squeeze(cat(1,temp{:}));
    temp=squeeze(mean(temp,2));
    


    [x,y,z]=ellipsoid_Hanlin(mean(temp),cov(temp),1,25);
    surf(x,y,z, 'FaceColor', 'none','EdgeColor',clist(i,:), 'EdgeAlpha', 0.7);
    hold on;
end
view([49 21])
axis tight 
% zlim([-1 3]);
zticks([-1 0 1 2 3])
% xlim([-3 5]);
xticks([-2 0 2 4])
% ylim([-3.5 4]);
yticks([-2 0 2 4])
set(gca,'FontSize',18)

xlabel('CP-1','Rotation',-15)
ylabel('CP-2','Rotation',20)
zlabel('CP-3')
thisA=gca;
axes('Position',[0.23    0.0500    0.7188    0.9238])
c=colorbar;
colorList=hsv(16);
colorList(5,:) = [121 130 56]/255;
colorList(7,:) = [1 121 111]/255;
colormap(colorList);
ylabel(c,'degree')
c.Ticks=[(1/16):(1/16):(1)]-1/32;
c.TickLabels=num2cell(reshape(reshape(0:22.5:(360-22.5),8,2)',[],1));
c.Direction='reverse';
c.Position=[0.82    0.03    0.0306    0.88]
c.FontSize=17;
set(thisA,'Position',[0.1041    0.1105    0.6803    0.8145])
axis off
saveas(gcf,'results/Fig4c-right.png')
% exportgraphics(gcf,'results/Fig4c-right.pdf','ContentType','vector')
close(gcf)

%% plot Fig4d-this neuron
figure
scatter(LDA_decodePower_trialS,LDA_15CV_Store,40,'ko','filled');axis equal;hold on;
axisMax=0.38;chance=1/16;scatter(chance,chance,40,'bo','filled')
axis([0 axisMax 0 axisMax]); plot(0:0.01:axisMax,0:0.01:axisMax,'b--')
set(gca,'FontSize',16)
xticks([0:0.1:0.4]);yticks([0:0.1:0.4])
xlabel('accuracy, rates ')
ylabel('accuracy, components')
legend('this neuron','chance','equality')
saveas(gcf,'results/Fig4d-One-Neuron.png')
close(gcf)

%% plot Fig4d
load('NI-Fig4d.mat')
figure
scatter(rateCodeDecode, temporalCodeDecode,'k.');axis equal;hold on;
[mean(rateCodeDecode) sem(rateCodeDecode)]
[mean(temporalCodeDecode) sem(temporalCodeDecode)]
axisMax=xlim;temp=ylim;axisMax=max(axisMax(2),temp(2));
chance = 1/100;
scatter(chance,chance,40,'bo','filled')
axis([0 axisMax 0 axisMax]); plot(0:0.01:axisMax,0:0.01:axisMax,'b--')
[h,p,ci_t2,stats_t2]=ttest(rateCodeDecode,temporalCodeDecode);
set(gca,'FontSize',16)
xtickarray=[0:0.02:axisMax];
xticks(xtickarray);yticks(xtickarray)
xlabel('accuracy, rates ')
ylabel('accuracy, components')
text(range(xlim)*0.5,range(xlim)*0.1,sprintf('p = %.2g',p),'VerticalAlignment','bottom','FontSize',14)
title('NI','FontWeight','normal')
saveas(gcf,'results/Fig4d.png')
close(gcf)

%% plot Fig4e,f,g
stimNames = {'DG','SG','RFG','NI'};
for nVS=1:4
load([stimNames{nVS} '-CO.Fig4.mat'])
figure('Position',[ 343    43   432   410])
scatter(atanh(staFR), atanh(staCP),'bo','MarkerFaceColor','b','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.3);axis equal;hold on;
axisMax=xlim;temp=ylim;axisMax=max(axisMax(2),temp(2));
cla
plot(0:0.01:axisMax,0:0.01:axisMax,'k--');hold on;
scatter(atanh(staFR), atanh(staCP),'bo','MarkerFaceColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.3);axis equal;hold on;

axis equal
box off
axis([0 axisMax 0 axisMax]); 
xticks([0.5:0.5:2.5]);yticks([0.5:0.5:2.5])
if nVS==1
    xticks([0.5:0.5:3.5]);yticks([0.5:0.5:3.5])
end
var1=atanh(staFR); 
var2=atanh(staCP);


[h,p,ci_t2,stats_t2]=ttest(var1(~isnan(var1) & ~isnan(var2)),var2(~isnan(var1) & ~isnan(var2)) );
set(gca,'FontSize',17)
xlabel('rates ')
ylabel('components')
text(range(xlim)*0.5,range(xlim)*0.1,sprintf('p = %.2g',p),'VerticalAlignment','bottom','FontSize',15)
title(stimNames{nVS},'FontWeight','normal')

saveas(gcf,['results/Fig4e-' stimNames{nVS}  '.png'])
close(gcf)

% plot Fig4f
h=figure('Position',[1144         221         300         405]);
global xaxisshaded
xaxisshaded=0:14;
hanlinErrorBar(driftRateMG,'b')
    hold on;
     hanlinErrorBar(driftRateCPMG,'r')
    xlabel('interval (days)')
    ylabel('representational drift index')
    xaxisshaded=[];
xlim([-0.5 14.5])
l=legend('rate','components','Color','none','EdgeColor','none','Location',[0.2616 0.8237 0.5222 0.1219]);
set(gca,'FontSize',16)
title([stimNames{nVS} '-all'],'FontWeight','normal')
yl=ylim;

switch nVS
    case 1
ylim([0 0.176])
   case 2
ylim([0 0.145])
case 3
ylim([0 0.148])
 case 4
ylim([0 0.2070])
end
yl=ylim;

var1=driftRateCPMG;
var2=driftRateMG;

valid=any(~isnan(var1),2) | any(~isnan(var2),2) ;

var1=var1(valid,:);
var2=var2(valid,:);

lme=hanlinLME([{var1} {var2}],0);



str1=sprintf('p(same slope)\n=%.2g',double(lme{4}.Coefficients(4,6)));

if double(lme{4}.Coefficients(4,6))<realmin
str1=sprintf('p(same slope)\n<%.2g',realmin);
end

text(0.25,yl(2)*0.81,str1,'FontSize',14,'Color','k')
saveas(gcf,['results/Fig4f-' stimNames{nVS}  '.png'])
close(gcf)

% Fig4g

relThres=75;% number between 50-100
partString={'-least-reliable-'};
part=1;
relPart=relvalue<prctile(relvalue,100-relThres);

h=figure('Position',[1144         221         300         405]);
global xaxisshaded
xaxisshaded=0:14;
hanlinErrorBar(driftRateMG(relPart,:),'b')
    hold on;
     hanlinErrorBar(driftRateCPMG(relPart,:),'r')
    xlabel('interval (days)')
    ylabel('representational drift index')
    xaxisshaded=[];
xlim([-0.5 14.5])
l=legend('rate','components','Color','none','EdgeColor','none','Location',[0.2738 0.8191 0.5222 0.1198]);
set(gca,'FontSize',16)
title([stimNames{nVS} partString{part} num2str(100-relThres) '%'],'FontWeight','normal')


yl=ylim;
switch nVS
    case 1
ylim([0 0.23])
   case 2
ylim([0 0.173])
case 3
ylim([0 0.268])
 case 4
ylim([0 0.203])
end
yl=ylim;
var1=driftRateCPMG(relPart,:);
var2=driftRateMG(relPart,:);

valid=any(~isnan(var1),2) | any(~isnan(var2),2) ;

var1=var1(valid,:);
var2=var2(valid,:);

lme=hanlinLME([{var1} {var2}],0);

str1=sprintf('p(same slope)\n=%.2g',double(lme{4}.Coefficients(4,6)));
if double(lme{4}.Coefficients(4,6))<realmin
str1=sprintf('p(same slope)\n<%.2g',realmin);
end
text(0.25,yl(2)*0.8,str1,'FontSize',14,'Color','k')
saveas(gcf,['results/Fig4g-' stimNames{nVS}  '.png'])
close(gcf)

end
