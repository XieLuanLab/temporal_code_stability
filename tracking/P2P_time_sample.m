function output = P2P_time_sample(waveform, negative)
% assume peak happens after valley 
waveform=resample(waveform,10,1);
if negative==-1
    minSample = find(waveform==min(waveform));
    minSample=minSample(1);
    maxSample = find(waveform(minSample:end)==max(waveform(minSample:end)))+minSample-1;
    thres=(max(waveform(minSample:end))-min(waveform))*0.95+min(waveform);
    output = sum(waveform(minSample:maxSample)<=thres);
else
    minSample = find(waveform==max(waveform));
    minSample=minSample(1);
    maxSample = find(waveform(minSample:end)==min(waveform(minSample:end)))+minSample-1;
    thres=max(waveform)-(max(waveform)-min(waveform(minSample:end)))*0.95;
    output = sum(waveform(minSample:maxSample)>=thres);
end
output=output/10;
end