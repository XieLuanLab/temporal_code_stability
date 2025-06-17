temp=[];temp.placeholder=1;
TrackSum{end+1,1}=temp;
switch nVS
    case 1
        temp=load(ustb_fileName,'VSFeat_Store','USTB_DG');
        TrackSum{end}.DG=cellfun(@(x) x(:,:,1:floor(VSTCB.DG(3)/2/10)),temp.USTB_DG,'UniformOutput',false);
        if timeAverage
            TrackSum{end}.DG=cellfun(@(x) squeeze(mean(double(x),3)),TrackSum{end}.DG,'UniformOutput',false);
        end
        if trialAverage
            TrackSum{end}.DG=cellfun(@(x) squeeze(mean(double(x),1)),TrackSum{end}.DG,'UniformOutput',false);
        end
    case 2
        temp=load(ustb_fileName,'VSFeat_Store','USTB_SG');
        TrackSum{end}.SG=cellfun(@(x) x(:,:,1:floor(VSTCB.SG(3)/2/10)),temp.USTB_SG,'UniformOutput',false);
        if timeAverage
            TrackSum{end}.SG=cellfun(@(x) squeeze(mean(double(x),3)),TrackSum{end}.SG,'UniformOutput',false);
        end
        if trialAverage
            TrackSum{end}.SG=cellfun(@(x) squeeze(mean(double(x),1)),TrackSum{end}.SG,'UniformOutput',false);
        end
    case 3
        temp=load(ustb_fileName,'VSFeat_Store','USTB_RFG');
        TrackSum{end}.RFG=cellfun(@(x) x(:,:,1:floor(VSTCB.RFG(3)/2/10)),temp.USTB_RFG,'UniformOutput',false);
        if timeAverage
            TrackSum{end}.RFG=cellfun(@(x) squeeze(mean(double(x),3)),TrackSum{end}.RFG,'UniformOutput',false);
        end
        if trialAverage
            TrackSum{end}.RFG=cellfun(@(x) squeeze(mean(double(x),1)),TrackSum{end}.RFG,'UniformOutput',false);
        end
    case 5
        temp=load(ustb_fileName,'VSFeat_Store','USTB_IJO');
        TrackSum{end}.IJO=cellfun(@(x) x(:,:,1:floor(VSTCB.IJO(3)/2/10)),temp.USTB_IJO,'UniformOutput',false);
        if timeAverage
            TrackSum{end}.IJO=cellfun(@(x) squeeze(mean(double(x),3)),TrackSum{end}.IJO,'UniformOutput',false);
        end
        if trialAverage
            TrackSum{end}.IJO=cellfun(@(x) squeeze(mean(double(x),1)),TrackSum{end}.IJO,'UniformOutput',false);
        end
end
