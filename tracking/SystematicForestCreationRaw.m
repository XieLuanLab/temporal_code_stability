%% generate a 4 dimension forest fully made of hierachrchical trees.
function [rawZ]= SystematicForestCreationRaw...
    (DistFeatForTreeBuild,XCorrDistFeatForTreeBuild,mask)
'generate a 4 dimension forest'
% input =  flattened waveform  Xnew
% output = list(cell array) of trees
distMatrixTreeA = cell(2,5,3,4); % so you have 4 dimensions to describe each tree 
rawZ = cell(1,numel(DistFeatForTreeBuild)+numel(XCorrDistFeatForTreeBuild));
% dimension 1  = [ raw waveform , first derivative ] 
% dimension 2 =  [distance type used ]
% dimension 3 =  [type of feature used , PCA, UMAP]
% dimension 4(5) = [type of link used ] % since dimension 4 is a singleton dimension. not
% used  here. 
treeNum=1;
for wav = 1:2 % dimension 1
    for interest = 1:5 % dimension 2
        for combine_dist = 1:3 % dimension 3 distCurrent  = distMatrix(:,:,combine_dist);
            Y =DistFeatForTreeBuild{wav,interest,combine_dist};
            Y(mask==0)=Y(mask==0)+max(Y(:)); % some hypothesis that the relative order is still important but incoporting prior is also important. 
%             Y(mask==1)=min(Y(:)); % some furter assumption that known clusters link to themselves first. 0620
            Y = Y - diag(diag(Y)); % make sure the diagnoal is zero
            yOut = squareform(Y,'tovector');
            % try hierarchical --------------------------------------------------------------------------
            Z = cell(4,1);
            Z{1} = linkage(yOut,'average');
            Z{2} = linkage(yOut,'single');
            Z{3} = linkage(yOut,'complete');
            Z{4} = linkage(yOut,'weighted');
            for linktype = 1:4 % dimension 4(5)
                rawZ{treeNum}=Z{linktype};
                treeNum=treeNum+1;
                generatingTree = [wav,interest,combine_dist,linktype];
                %    distMatrixTree{wav,interest,combine_dist,inter,linktype} = Z{linktype};
                %    distMatrixIdx{wav,interest,combine_dist,inter,linktype}=[wav interest combine_dist inter linktype];
            end
        end
    end
end


for i=1:2
    for j =1:6
Y =XCorrDistFeatForTreeBuild{i,j};
            Y(mask==0)=Y(mask==0)+max(Y(:)); % some hypothesis that the relative order is still important but incoporting prior is also important. 
            Y = Y - diag(diag(Y)); % make sure the diagnoal is zero
            yOut = squareform(Y,'tovector');
            % try hierarchical --------------------------------------------------------------------------
            Z = cell(4,1);
            Z{1} = linkage(yOut,'average');
            Z{2} = linkage(yOut,'single');
            Z{3} = linkage(yOut,'complete');
            Z{4} = linkage(yOut,'weighted');
            for linktype = 1:4 % dimension 4(5)
                rawZ{treeNum}=Z{linktype};
                treeNum=treeNum+1;
                generatingTree = [i j];
                %    distMatrixTree{wav,interest,combine_dist,inter,linktype} = Z{linktype};
                %    distMatrixIdx{wav,interest,combine_dist,inter,linktype}=[wav interest combine_dist inter linktype];
            end    
    end
end
% rawZIdx = reshape(distMatrixIdx,1,numel(distMatrixIdx));
'done'
