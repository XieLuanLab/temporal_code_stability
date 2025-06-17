
for nVS=[1 2 3 5]
switch nVS
    case 2
        nPattern=30;
    case 1
        nPattern=16;
    case 5
        nPattern=100;
    case 3
        nPattern=81;
end
basePath=pwd;
addpath(genpath(basePath))


stimNames = {'DG','SG','RFG','','IJO'};
greaterThanNd=15;
zmode=0;
sampleEveryNTrial=2; % so you have even and odd trials
umapThres=85;


nthComp=1:3;
perDayCVcomp=false;
ZeroMean=0;
DivideSum=1;
kernalWidth=24;
kfold=5;
updatedRate=1;
prepareUMAP_FRdata()
updatedComponents=1;
prepareUMAP_CPdata();
%
samePatternAcrossDays=cell(5,15);
diffPatternWithinDays=cell(5,15);
diffPatternWithinDaysNN=cell(5,15,3);% give me nearest 1, 5, 10.
PatternIDperItem=repmat(1:nPattern,2*15,1)';
PatternIDperItem=PatternIDperItem(:);
PatternIDperItemDist=pdist2(PatternIDperItem,PatternIDperItem);
PatternIDperItemDistList=pdist(PatternIDperItem);
PatternIDperItemList=PatternIDperItem+zeros(size(PatternIDperItem))';
PatternIDperItemList=return_ltri(PatternIDperItemList')';
PatternIDperItemList2=PatternIDperItem+zeros(size(PatternIDperItem))';
PatternIDperItemList2=return_ltri(PatternIDperItemList2)';

dayNumberPerItem=repmat(1:15,nPattern*2,1);
dayNumberPerItem=dayNumberPerItem(:);
dayNumberPerItemDist=pdist2(dayNumberPerItem,dayNumberPerItem);
dayNumberPerItemDistList=pdist(dayNumberPerItem);

lastDayNum=max(dayNumberPerItem,dayNumberPerItem'); % added lines 
lastDayNum=lastDayNum-diag(diag(lastDayNum));
lastDayNumPerItemDistList=squareform(lastDayNum,'tovector');

firstDayNum=min(dayNumberPerItem,dayNumberPerItem'); % added lines 
firstDayNum=firstDayNum-diag(diag(firstDayNum));
firstDayNumPerItemDistList=squareform(firstDayNum,'tovector');

dayNumberPerItemList=dayNumberPerItem+zeros(size(dayNumberPerItem))';
dayNumberPerItemList=return_ltri(dayNumberPerItemList')';


typeStrList={'FR','CP'};
for tt=1:2
typeStr=typeStrList{tt};%input('FR or CP','s');
for ani = 1:5

data2 = cat(1,designMatrixCP{ani,:});
 
data1=cat(1,designMatrix{ani,:});    

if  strcmp(typeStr,'FR')
        data=data1;
   elseif strcmp(typeStr,'CP')
       data=data2;
   end

switch nVS
    case 2
temp=pdist2(data,data,@distfunEUCL_nan);%SG-euclidean=-134-0.7
temp=(temp+temp')/2;
temp=temp-diag(diag(temp));
[scored3, umap, ~, ~]=run_umap(temp,'randomize',false,'n_neighbors',134,'n_components',2,'min_dist',0.7,'metric','precomputed','cluster_output','none','verbose','none');%sg
    case 1
temp=1-corr(data','rows','pairwise');%DG-corr-149-0.9
temp=(temp+temp')/2;
temp=temp-diag(diag(temp));
[scored3, umap, ~,~]=run_umap(temp,'randomize',false,'n_neighbors',149,'n_components',2,'min_dist',0.9,'metric','precomputed','cluster_output','none','verbose','none');%DG
    case 5
temp=1-corr(data','rows','pairwise');%NI-corr-119-0.9
temp=(temp+temp')/2;
temp=temp-diag(diag(temp));
[scored3, umap, ~, ~]=run_umap(temp,'randomize',false,'n_neighbors',119,'n_components',2,'min_dist',0.9,'metric','precomputed','cluster_output','none','verbose','none');%NI
    case 3
data=pdist2(data,data,@distfunCosNaN);%RFG-cos-134-0.9
data=data-eye(size(data)).*data; data=(data+data')/2;
[scored3, umap, ~, ~]=run_umap(data,'randomize',false,'n_neighbors',134,'n_components',2,'min_dist',0.9,'metric','precomputed','cluster_output','none','verbose','none'); %rfG
end

patternDistMatrixList=pdist(scored3);
daysAcc=30;
    temp=1;
    for d=1:15
        targetDist=patternDistMatrixList(dayNumberPerItemDistList==0 & PatternIDperItemDistList~=0 & dayNumberPerItemList==d);
        targetDistPatternID=PatternIDperItemList(dayNumberPerItemDistList==0 & PatternIDperItemDistList~=0 & dayNumberPerItemList==d);
        targetDistPatternID2=PatternIDperItemList2(dayNumberPerItemDistList==0 & PatternIDperItemDistList~=0 & dayNumberPerItemList==d);
        Pat2PatMat=zeros(nPattern);
        for item=1:numel(targetDist)
        Pat2PatMat(targetDistPatternID(item),targetDistPatternID2(item))=Pat2PatMat(targetDistPatternID(item),targetDistPatternID2(item))+targetDist(item);
        end
        Pat2PatMat=Pat2PatMat+Pat2PatMat';
        Pat2PatMat=Pat2PatMat+max(Pat2PatMat(:))*eye(nPattern);
        target=mink(Pat2PatMat,10);
        diffPatternWithinDaysNN{ani,d,3}=target(:)/temp/4;
        target=target(1:5,:);
        diffPatternWithinDaysNN{ani,d,2}=target(:)/temp/4;
        target=target(1:1,:);
        diffPatternWithinDaysNN{ani,d,1}=target(:)/temp/4;
        
        diffPatternWithinDays{ani,d}=targetDist'/temp;
    end


temptemp=cat(1,diffPatternWithinDaysNN{ani,:,2}); % target NN=5
temptemp=squeeze(mean(reshape(temptemp,5,nPattern,[]),1));
    for d=0:14
        targetDist=patternDistMatrixList(dayNumberPerItemDistList==d & PatternIDperItemDistList==0);
        targetDistPatternID=PatternIDperItemList(dayNumberPerItemDistList==d & PatternIDperItemDistList==0);
        targetDistlastDayID=lastDayNumPerItemDistList(dayNumberPerItemDistList==d & PatternIDperItemDistList==0);
        targetDistfirstDayID=firstDayNumPerItemDistList(dayNumberPerItemDistList==d & PatternIDperItemDistList==0);
        samePatternAcrossDays{ani,d+1}=arrayfun(@(x) mean(targetDist(targetDistPatternID==x)...
            ./(temptemp(x,targetDistlastDayID(targetDistPatternID==x)) +temptemp(x,targetDistfirstDayID(targetDistPatternID==x)))*2  ), 1:nPattern)';
    % equation 8 
    end

if ani==3 && nVS==5
h=figure();
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 13.2;
h.Position(4)         = 12;
h.PaperPositionMode   = 'auto';
colorList=distinguishable_colors(100,[1 1 1]);
hold on;
for img=1:100
plot(scored3(img:100:end,1),scored3(img:100:end,2),'Color',colorList(img,:),'LineWidth',0.5)
end
scatter(scored3(:,1),scored3(:,2),6,repmat( colorList(1:100,:),daysAcc,1),'o','filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
xlabel('umap-1');%xticks([-10 0 10])
ylabel('umap-2');%xticks([-10 0 10])
axis tight
set(gca,'FontSize',12)
c=colorbar;
colormap( colorList)
ylabel(c,'Image ID','FontSize',10)
c.Ticks=[0.1:0.1:1.0]-1/200;
c.TickLabels=num2cellstr([10:10:100]);
title([ 'NI-' typeStr],'FontWeight','normal','FontSize',12)
exportgraphics(h,[stimNames{nVS} '-n-' typeStr '.pdf'],'ContentType','vector')
pause(1);close(h)
end


if ani==1 && nVS==3
h=figure('Position',[100        100        796   680]) ;
hold on
for i=1:9*daysAcc
 patch([-scored3(9*(i-1)+1:i*9,2); nan],[-scored3(9*(i-1)+1:i*9,1); nan],mod(i-1,9)*1.5+1*[1:9 nan],'EdgeColor',[0.8 0.8 0.8],'FaceColor','none','lineWidth',.5,'EdgeAlpha',0.2)
end
cmap=repmat(turbo(9),9,1);
for pat=1:81
     temp = scored3(pat:81:end,:);
l=plot(-temp(:,2),-temp(:,1),'Color',cmap(pat,:));
l.Color(4)=0.6;
end
hold on
scatter(-scored3(:,2),-scored3(:,1),15,repmat(turbo(9),9*daysAcc,1),'o','filled','MarkerFaceAlpha',0.6)
set(gca,'FontSize',16)
xlabel('umap-2')
ylabel('umap-1')
colormap((turbo(9)))
c=colorbar;
title('')
clim([0 0.9])
c.Ticks=[0.05:0.2:0.85];c.Ticks=[0.05:0.1:0.85];
a=7.4;
[x,y]=meshgrid(-4:4,-4:4);
x=x*3.7;
y=flip(y)*3.7;
d=20.5;
azim = 45-atand(x/d);
elev = atand(y./sqrt(x.^2+d^2));
rfsd = (0.5*a./sqrt(x.^2+y.^2+d^2))/pi*180;
c.TickLabels=num2cellstr(([azim(1,1:2:end)]));c.TickLabels=num2cellstr(([azim(1,1:1:end)]));
c.Location='north';
ylabel(c,'Azim (degree)','FontSize',14);

if strcmp(typeStr,'CP')
for i=1:9
    text(i*1.5-5,-8+1.5*i,sprintf('Elev %0.2f',elev(10-i,1)),'FontSize',16)
end
end
title([stimNames{nVS} '-' typeStr],'FontWeight','normal','FontSize',16)
exportgraphics(h,[stimNames{nVS} '-n-' typeStr '.png'])
pause(1);close(h);
end

if ani==1 && nVS==2
h=figure('Position',[100        100        796         680]);
baseMap=hsv(7);  baseMap=baseMap([1 2 3 5 6 7],:);
  cmap=repmat(baseMap,5,1);
  hold on;
  for pat=1:30
     temp = scored3(pat:30:end,:);
l=plot(temp(:,1),temp(:,2),'Color',cmap(pat,:));
l.Color(4)=0.7;
  end
   markers={'o','v','p','^','diamond'};
allMaps=repmat(baseMap,5*daysAcc,1);
for fre=1:5
    ms=20;
    if fre==3
    ms=40;
    end
    for d = 1:daysAcc
      valid=  sri(30,d,(fre-1)*6+1:fre*6);
score3 =spline(0:length(scored3(valid,1))-1,scored3(valid,:)',0:0.1:length(scored3(valid,1))-0.1)';
patch([score3(1:end-10,1); nan],[score3(1:end-10,2); nan],19*ones(size([score3(1:end-10,1);nan])),'EdgeColor',[0.8 0.8 0.8],'FaceColor','none','lineWidth',.5,'EdgeAlpha',0.2)
scatter(scored3(valid,1),scored3(valid,2),ms,allMaps(valid,:),markers{fre},'LineWidth',1)
hold on
    end
end
colormap(baseMap);
c=colorbar;
ylabel(c,'degree')
c.Ticks=([(1/6):(1/6):(1)]-1/12)*2+18;
c.TickLabels=num2cellstr(0:30:(180-30));
c.Direction='reverse';
set(gca,'FontSize',16)
xlabel('umap-1')
ylabel('umap-2')
title([stimNames{nVS} '-' typeStr],'FontWeight','normal','FontSize',16)
exportgraphics(h,[stimNames{nVS} '-n-' typeStr '.png'])
pause(1);close(h);
end
if ani==3 && nVS==1
h=figure();
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 13.2;
h.Position(4)         = 12;
h.PaperPositionMode   = 'auto';
colorList=hsv(16);
colorList=[colorList(1:2:end,:); colorList(2:2:end,:)];
colorList(3,:) = [121 130 56]/255;
colorList(4,:) = [1 121 111]/255;
score3 =spline(0:length(scored3)-1,scored3',0:0.1:length(scored3)-0.1)';
 patch([score3(10:end-9,1); nan],[score3(10:end-9,2); nan],19*ones(size([score3(10:end-9,2);nan])),'EdgeColor',[0.8 0.8 0.8],'FaceColor','none','lineWidth',.5,'EdgeAlpha',0.2)
hold on;
for pat=1:16
     temp = scored3(pat:16:end,:);
     l=plot(temp(:,1),temp(:,2),'Color',colorList(pat,:),'lineWidth',0.5);
     l.Color(4)=0.4;
end
 scatter(scored3(:,1),scored3(:,2),9,repmat(colorList,daysAcc,1),'o','filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
axis tight
xlabel('umap-1');%xticks([-10 0 10])
ylabel('umap-2');%yticks([-4 0 4])
set(gca,'FontSize',12)
title([stimNames{nVS} '-' typeStr],'FontWeight','normal','FontSize',12)
exportgraphics(h,[stimNames{nVS} '-n-' typeStr '.pdf'],'ContentType','vector')
pause(1);close(h);
end

end
tlstr=[basePath filesep stimNames{nVS} '-n-umap-' typeStr];
save(tlstr,'samePatternAcrossDays','diffPatternWithinDays',"diffPatternWithinDaysNN");
end
end
%% OOO actual statictis
for nVS=[1 2 3 5]

typeStr='CP';
tlstr=[basePath filesep stimNames{nVS} '-n-umap-' typeStr];

load(tlstr)
driftRateCP=arrayfun(@(x) cat(1,samePatternAcrossDays{:,x}),1:15,'UniformOutput',false);
driftRateCP=cell2mat(driftRateCP);

typeStr='FR';
tlstr=[basePath filesep stimNames{nVS} '-n-umap-' typeStr];
load(tlstr)
driftRate=arrayfun(@(x) cat(1,samePatternAcrossDays{:,x}),1:15,'UniformOutput',false);
driftRate=cell2mat(driftRate);


h=figure;
h.Units               = 'centimeters';
h.Position(1)         = 3;
h.Position(2)         = 3;
h.Position(3)         = 3.8*2*1.05;
h.Position(4)         = 4.65*2;
h.PaperPositionMode   = 'auto';

global xaxisshaded
xaxisshaded=0:14;
plot([-0.5 14.5],[1 1],'k--')
hanlinErrorBar(driftRate,'b')
    hold on;
     hanlinErrorBar(driftRateCP,'r')
    xlabel('interval (days)')
   
ylabel({'pattern drift against','nearest 5 patterns'})
xlim([-0.5 14.5])

switch nVS
    case 1
ylim([0.16 1.36])
   case 2
ylim([0.18 2.25])
case 3
ylim([0.5 7])
 case 5
ylim([0.25 2.4])
end
yl=ylim;

l=legend('','rate','components','Color','none','EdgeColor','none','Location',[0.275 0.8237 0.5222 0.1219]);

title(stimNames{nVS},'FontWeight','normal')
if nVS==5
    title('NI','FontWeight','normal')
end
set(gca,'FontSize',15)

yl=ylim;

var1=driftRateCP;
var2=driftRate;
[r,c]=find(isnan(var1)~=isnan(var2));
valid=any(~isnan(var1),2) | any(~isnan(var2),2) ;

var1=var1(valid,:);
var2=var2(valid,:);
[nanmean(var1(:,end)) sem(var1(:,end))];
[nanmean(var2(:,end)) sem(var2(:,end))];

lme=hanlinLME([{var1} {var2}],0);
str1=sprintf('p(same slope)\n=%.2g',double(lme{4}.Coefficients(4,6)));

if double(lme{4}.Coefficients(4,6))<realmin
str1=sprintf('p(same slope)\n<%.2g',realmin);
end
box off

set(findall(h,'-property','FontSize'),'FontSize',14)
text(0.25,range(yl)*0.81+yl(1),str1,'FontSize',12,'Color','k')

set(findall(h,'Type','Legend','-property','FontSize'),'FontSize',12)

exportgraphics(h,[stimNames{nVS} '-umapSN5.pdf'],'ContentType','vector')
pause(1);close(h);
end
