% clear all; close all; % delete(gcp);
standardFs = 25;
files1=dir('*day*.mat');

if skipOneFile==1
files = files1(1:2:end); % if 
else
files = files1; 
end

if skipOneFile==1 && ~strcmp(files(end).name(13:14),'23') && ~strcmp(files(end).name,files1(end).name)
    files = files1([1:2:end end]); % if 
end
%% importData into cells.
AllData = cell(size(files,1),1);
for i=1:size(AllData,1)
AllData{i}=importdata(files(i).name);
AllData{i}.ValleyAcrossTime=[];
AllData{i}.PeakAcrossTime=[];    
end
MaxNumelChs =min(cellfun(@(x) size(x.Dat_V_Map,1), AllData));

dayNumForEachUnit  = cell(size(AllData,1),1);
for i = 1:size(AllData,1)
dayNumForEachUnit{i} = i*ones(numel(AllData{i}.waveform),1);
end
dayNumForEachUnit = vertcat(dayNumForEachUnit{:});


whichChsFeature = cell(size(AllData,1),1);
for i=1:size(AllData,1)
whichChsFeature{i}=repmat({AllData{i}.Dat_V_Map(:,2)},numel(AllData{i}.waveform),1);
end
whichChsFeature = cat(2,whichChsFeature(:));
whichChsFeature = vertcat(whichChsFeature{:});

ChMapFeature = cell(size(AllData,1),1);
for i=1:size(AllData,1)
ChMapFeature{i}=AllData{i}.Ch_Map;
end

DataForFeature = AllData{1}.waveform';
for i=2:size(AllData,1)
DataForFeature=[DataForFeature;(AllData{i}.waveform)'];
end


ValleyList = min(DataForFeature{1}');
minCh = find(ValleyList==min(ValleyList));
for i=2:size(DataForFeature,1)
ValleyList = min(DataForFeature{i}');
listOfmin = find(ValleyList==min(ValleyList));
minCh = [minCh;listOfmin(1)];
end



FRForFeature = AllData{1}.Avg_FR';
for i=2:size(AllData,1)
FRForFeature=[FRForFeature;AllData{i}.Avg_FR'];
end


StdForFeature = AllData{1}.waveformStd';
for i=2:size(AllData,1)
StdForFeature=[StdForFeature;AllData{i}.waveformStd'];
end


%% 


P2PForFeature = AllData{1}.P2P';
for i=2:size(AllData,1)
P2PForFeature=[P2PForFeature;AllData{i}.P2P'];
end


PeakChList = zeros(numel(P2PForFeature),1);
PosNegUnit = zeros(numel(P2PForFeature),1);

for i=1:numel(P2PForFeature)
pl = P2PForFeature{i};
[~,I1]=sort(pl,'desc');
PeakChList(i)=I1(1);
waveform=DataForFeature{i}(PeakChList(i),:);
if abs(min(waveform))< max(waveform)
PosNegUnit(i)=1;
else
PosNegUnit(i)=-1;

end
end

FWHM_summary = zeros(numel(P2PForFeature),1);
valley2pkTime_summary = zeros(numel(P2PForFeature),1);
pkCh_pkTime_summary = zeros(numel(P2PForFeature),1);
for i=1:numel(P2PForFeature)
FWHM_summary(i) = FWHM(double(DataForFeature{i}(PeakChList(i),:)),PosNegUnit(i));
valley2pkTime_summary(i) = P2P_time_sample(double(DataForFeature{i}(PeakChList(i),:)),PosNegUnit(i));
[~,pkCh_pkTime_summary(i)]= find(abs(DataForFeature{i})==max(abs(DataForFeature{i}),[],'all'),1);
end



top3CV=[];
top3MeanMinusStd=[];

for i=1:numel(DataForFeature)
    i
[~, I] = sort(P2PForFeature{i},'desc');
for top = 1:min(3,MaxNumelChs)
if abs(min(DataForFeature{i}(I(top),:)))>max(DataForFeature{i}(I(top),:))
[minV, minI] = min(DataForFeature{i}(I(top),:));
minI=minI(1);
[minstd] = min(StdForFeature{i}(I(top),minI)); 
top3CV(i,top) = minstd/abs(minV);
top3MeanMinusStd(i,top) = abs(minV)-minstd;
else

    [minV, minI] = max(DataForFeature{i}(I(top),:));
minI=minI(1);
[minstd] = max(StdForFeature{i}(I(top),minI)); 
top3CV(i,top) = minstd/abs(minV);
top3MeanMinusStd(i,top) = abs(minV)-minstd;
end

end
end


clear p2pDist p2pFit crossP
%parpool('local', 16)
p2pDist = zeros(numel(DataForFeature),MaxNumelChs);
crossP  = cell(numel(DataForFeature),1);
for i=1:numel(DataForFeature)
%     i
[AmpList, I] = sort(P2PForFeature{i},'desc');
top = 1;
if abs(min(DataForFeature{i}(I(top),:)))>max(DataForFeature{i}(I(top),:))
[minV, minI] = min(DataForFeature{i}(I(top),:));
p2pDist(i,top)= minV;
crossP{i}= min(DataForFeature{i},[],2);
else
[minV, minI] = max(DataForFeature{i}(I(top),:));
p2pDist(i,top)= minV;
crossP{i}= max(DataForFeature{i},[],2);
end

for top = 2:min(size(DataForFeature{i},1),20)
if p2pDist(i,1)<0
[minV, minI] = min(DataForFeature{i}(I(top),:));
p2pDist(i,top)= minV;
else
[minV, minI] = max(DataForFeature{i}(I(top),:));
p2pDist(i,top)= minV;  
end
end
end

timeForFeature = [];
Fs=AllData{1}.Fs
CCG_PerSession=cell(size(AllData,1),1); 
for i=1:size(AllData,1)
thisSesTimeForFeature = cellfun(@(x) double(x')/Fs, AllData{i}.time, 'UniformOutput', false);

timeForFeature=[timeForFeature thisSesTimeForFeature];
end

spikes.times=timeForFeature
spikes.numcells = numel(spikes.times);
spikes.total = cellfun(@(x) numel(x),spikes.times);
spikes.shankID = ones(spikes.numcells,1);
spikes.cluID = [1:spikes.numcells]';
acg_metrics = calc_ACG_metrics(spikes,Fs,'showFigures',false);
spikes=[];
timeForFeature=[];
normalized_p2pDist = abs(p2pDist(:,1:min(MaxNumelChs,20)))./abs(repmat(p2pDist(:,1),1,min(MaxNumelChs,20)));
sumAmp = mean(normalized_p2pDist(:,2:min(MaxNumelChs,20))'); 
% if lots of chs are turned off, this metrics could back fire, but I argue
% that it's more of electrode's problem, not of the metric itself
%% Noise Rejection
thres = [];
thres(:,1) = (abs(p2pDist(:,1))>=25); % 
thres(:,2) = FRForFeature>=0.01; % 
thres(:,3) = ~(top3CV(:,1)>2/3 & top3CV(:,2)>2/3 & top3MeanMinusStd(:,1)<40); % 
thres(:,4) = valley2pkTime_summary>=2.2*trueFs/standardFs & valley2pkTime_summary<=25.5*trueFs/standardFs; % 
thres(:,5) = FWHM_summary >=1.6*trueFs/standardFs &  FWHM_summary <=22.5*trueFs/standardFs;  % 
thres(:,6) = sumAmp<=2^(-0.05*min(MaxNumelChs,20)); % 
thres(:,7) = abs(pkCh_pkTime_summary-31)<=3*trueFs/standardFs; % 

remain_index = prod(thres',1); %
sum(remain_index)
sum(remain_index)/numel(DataForFeature)
result = AllData{end};

%% WithinSesCCG
CCG_PerSession=cell(size(AllData,1),1); 
WithinSesValidCellId = [];
for i=1:size(AllData,1)
thisSesTimeForFeature = cellfun(@(x) double(x')/Fs, AllData{i}.time, 'UniformOutput', false);
spikes=[]
thisSesRemain = remain_index(dayNumForEachUnit==i);
WithinSesValidCellId=[WithinSesValidCellId 1:sum(thisSesRemain)];
spikes.times=thisSesTimeForFeature(thisSesRemain==1);
spikes.numcells = numel(spikes.times);
spikes.total = cellfun(@(x) numel(x),spikes.times);
spikes.shankID = ones(spikes.numcells,1);
spikes.cluID = [1:spikes.numcells]';
mono_res = ce_MonoSynConvClick(spikes,'includeInhibitoryConnections',false); % detects the monosynaptic connections
CCG_PerSession{i}=single(mono_res.ccgR);
end
mono_res=[];
thisSesTimeForFeature=[];
spikes=[];
%% clean up 
AllData=[];
% 
%% what if you ch number differ across session, must find intersected chs across all 
if any(bsxfun(@minus,cat(3,ChMapFeature{:}),mean(cat(3,ChMapFeature{:}),3))~=0,'all')
    'alert'
finalMapAccepted = min(cat(3,ChMapFeature{:}),[],3);
[a,b]=sort(finalMapAccepted(:));
result.Ch_Map=finalMapAccepted;
new = zeros(4,8); % 
new(b)=1:32;
numZeros = sum(finalMapAccepted==0,'all');
new = new -numZeros;
new = new.*(new>0);
result.Ch_Map_2=new;

for i = 1:numel(DataForFeature)
    i
    extraChs = setdiff(whichChsFeature{i},finalMapAccepted);
    if ~isempty(extraChs)
        [~,~,ic]=intersect(extraChs,whichChsFeature{i});

    DataForFeature{i}(ic,:)=[];
    StdForFeature{i}(ic,:)=[];
    P2PForFeature{i}(ic)=[];
    crossP{i}(ic)=[];
    end
    
end
crossP = cat(2,crossP{:})';

PeakChList = zeros(numel(P2PForFeature),1);
PosNegUnit = zeros(numel(P2PForFeature),1);

for i=1:numel(P2PForFeature)
pl = P2PForFeature{i};
[~,I1]=sort(pl,'desc');
PeakChList(i)=I1(1);
waveform=DataForFeature{i}(PeakChList(i),:);
if abs(min(waveform))< max(waveform)
PosNegUnit(i)=1;
else
PosNegUnit(i)=-1;

end
end

FWHM_summary = zeros(numel(P2PForFeature),1);
valley2pkTime_summary = zeros(numel(P2PForFeature),1);
pkCh_pkTime_summary = zeros(numel(P2PForFeature),1);
for i=1:numel(P2PForFeature)
FWHM_summary(i) = FWHM(double(DataForFeature{i}(PeakChList(i),:)),PosNegUnit(i));
valley2pkTime_summary(i) = P2P_time_sample(double(DataForFeature{i}(PeakChList(i),:)),PosNegUnit(i));
[~,pkCh_pkTime_summary(i)]= find(abs(DataForFeature{i})==max(abs(DataForFeature{i}),[],'all'),1);
end


else
    crossP = cat(2,crossP{:})';
end
