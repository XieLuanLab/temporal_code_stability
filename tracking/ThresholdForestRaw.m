function [rawZ] = ThresholdForestRaw(DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,AllPartition)
%% generate a 4 dimension distanceMatrix , get 2sigma thresholds 
% input =  flattened waveform  Xnew
% output = FeatureThres
rawZ =[];
% dimension 1  = [ raw waveform , first derivative ] 
% dimension 2 =  [distance type used ]
% dimension 3 =  [type of feature used , PCA, UMAP]
% dimension 4(5) = [type of link used ] % since dimension 4 is a singleton dimension. not
tree = 0;
for wav = 1:2 % dimension 1 
for interest = 1:5 % dimension 2 
for combine_dist = 1:3 % dimension 3
distCurrent  = DistFeatForTreeBuild{wav,interest,combine_dist}; 
[Z] = Link4From4PartitionsDistanceMatrix(distCurrent,AllPartition(tree+1:tree+4));
rawZ = [rawZ;Z];
tree = tree+4
end
end
end

for i=1:2
    for j =1:6

distCurrent  = XCorrDistFeatForTreeBuild{i,j}; % this dist not regulated by MaskTrack_FarLocPerTopChXCorr because it has already by implicitly incoporatexdin calculation steps above.

[Z] = Link4From4PartitionsDistanceMatrix(distCurrent,AllPartition(tree+1:tree+4));
rawZ = [rawZ;Z];
tree = tree+4
    end
end

end


