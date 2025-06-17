function output = FWHM(waveform, negative)
waveform=resample(waveform,10,1);
if negative==-1
    minV = min(waveform);
    output = numel(waveform(waveform<=0.5*minV));
else
    minV = max(waveform);
    output = numel(waveform(waveform>=0.5*minV));
end
output=output/10;
end