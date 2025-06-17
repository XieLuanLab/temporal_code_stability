function h=plot_two_units(subplotRow,subplotCol,options,DataForFeatureAvgPer32days)
%              options.colorList=jet(numel(unique(dayNumForEachUnit))); 
%              options.Ch_Map=store1.Ch_Map-15;
%              options.Yoff= 250;
%              options.Xoff= 50;
%              options.Visible = 'on';
%              options.subplotRow = 3;
%              options.subplotCol = 2;
%              options.subplotContent = cell(options.subplotRow,options.subplotCol); % specify subplot content 
%              options.subplotContent{1,1}=leftChildComp;
%              options.subplotContent{1,2}=rightChildComp;
%              options.subplotContent{2,1}=[leftChildComp(1) rightChildComp(1)];
%              options.subplotContent{2,2}=[leftChildComp(end) rightChildComp(end)];
%              options.subplotContent{3,1}=[leftChildComp(1) rightChildComp(end)];
%              options.subplotContent{3,2}=[leftChildComp(end) rightChildComp(1)];


    list = options.subplotContent{subplotRow,subplotCol};
    Ch_Map = options.Ch_Map;
    try
    xCenter = options.xCenterTransferFunction(round(mean(options.xCenter(options.currentComposition))));
    Ch_Map = Ch_Map(:,xCenter-4:xCenter+4);
    catch
    end
    try
    alpha = options.alpha;
    catch
        alpha = 1;
    end
    

    colorList = options.colorList;

    
    Yoff = options.Yoff;
    Xoff = options.Xoff;
    
    numGx = size(Ch_Map,2);
    numGy = size(Ch_Map,1);
    t = 1:size(DataForFeatureAvgPer32days{1},2);
    
    elementPerTrace = numel(t);
    
    giantX = nan(1,numGy*numGx*numel(list)*(elementPerTrace+1));
    giantY = nan(1,numGy*numGx*numel(list)*(elementPerTrace+1)); % max amount of elements
    
    eidx = 0;
    for row=1:numGy 
        for col=1:numGx
                  
            ch = Ch_Map(row,col);
            if ch>0 % valid channel 
            for L = 1:numel(list)
                
                unit = list(L);
                    giantX(eidx+1:eidx+elementPerTrace)=t+Xoff*col;
                    giantY(eidx+1:eidx+elementPerTrace)=DataForFeatureAvgPer32days{unit}(ch,:)+Yoff*(numGy+1-row);
                    eidx = eidx + elementPerTrace +1 ; 
                    
%                     if alpha==1
%                     plot(t+Xoff*col,DataForFeatureAvgPer32days{unit}(ch,:)+Yoff*(numGy+1-row),'Color',colorList(1,:),'LineWidth',1);
%                     else
%                         
%                     plot1 = plot(t+Xoff*col,DataForFeatureAvgPer32days{unit}(ch,:)+Yoff*(numGy+1-row),'Color',colorList(1,:),'LineWidth',1);
%                     plot1.Color(4)=alpha;
%                     end 
%                 hold on;
                
            end
            end
        end 

    end
plot(giantX,giantY,'Color',colorList(1,:),'LineWidth',1);
hold on;
ylim([0 Yoff*numGy+100])
xlim([0 Xoff*numGx+50])
if subplotRow==1 && subplotCol==2
scatter(Xoff+options.dayNumForEachUnit(list)*Xoff/2,(1.5*Yoff+15)*ones(size(list)),90*ones(size(list)),'LineWidth',1.25,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.4,'Marker','square')
hold on;
list = setdiff(options.currentComposition,list);
scatter(Xoff+options.dayNumForEachUnit(list)*Xoff/2,(1.5*Yoff-15)*ones(size(list)),90*ones(size(list)),'LineWidth',1.25,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.4,'Marker','square')
end
% try
%     xCenter = options.xCenterTransferFunction(round(mean(options.xCenter(options.currentComposition))));
%     xlim([(xCenter-5)*Xoff (xCenter+5)*Xoff])
% catch
% end
% end
% set(gca,'Position',[ 0.0895    0.0892    0.8112    0.8612])
set(gca,'XTickLabel','')
axis off
% set(gca,'FontSize',14)
% xlabel('ms')
% ylabel('\muV')
% mkdir(['OnlyFour'])
% cd(['OnlyFour'])
% mkdir(['ROI16-ch-Wave'])
% cd(['ROI16-ch-Wave'])
% saveas(h,['RawWave-ch' num2str(chnls(i)) 'id-' num2str(img) '.png'])
% close(h)
% cd ..
end

% cd ..

