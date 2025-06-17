function plot_AcrossDay_units(DataForFeatureAvgPer32days,acg_metrics,dayNumForEachUnit,FRForFeature,options,shankStr)

    Ch_Map = options.Ch_Map;
%     try
%     xCenter = options.xCenterTransferFunction(round(mean(options.xCenter(options.currentComposition))));
%     Ch_Map = Ch_Map(:,xCenter-4:xCenter+4);
%     catch
%     end
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
    
   
parfor unit = 1:numel(DataForFeatureAvgPer32days) 
    h=figure('Position',[304         144        1500         900],'Visible','off');
    
for row=1:numGy 
        for col=1:numGx     
            ch = Ch_Map(row,col);
            if ch>0 % valid channel 
                    if alpha==1
                    plot(t+Xoff*col,DataForFeatureAvgPer32days{unit}(ch,:)+Yoff*(numGy+2-row),'Color','b','LineWidth',2.5);
                    else
                    plot1 = plot(t+Xoff*col,DataForFeatureAvgPer32days{unit}(ch,:)+Yoff*(numGy+2-row),'Color','b','LineWidth',2.5);
                    plot1.Color(4)=alpha;
                    end 
                hold on;              
            end
        end 
end
if max(DataForFeatureAvgPer32days{unit},[],'all')>200 || min(DataForFeatureAvgPer32days{unit},[],'all')<-400
for row=1:numGy 
        for col=1:numGx     
            ch = Ch_Map(row,col);
            if ch>0 % valid channel 
                    if alpha==1
                    plot(t+Xoff*col,DataForFeatureAvgPer32days{unit}(ch,:)/3+Yoff*(numGy+2-row),'Color', [.5 0 .5],'LineWidth',2.5);
                    else
                    plot1 = plot(t+Xoff*col,DataForFeatureAvgPer32days{unit}(ch,:)/3+Yoff*(numGy+2-row),'Color', [.5 0 .5],'LineWidth',2.5);
                    plot1.Color(4)=alpha;
                    end 
                hold on;              
            end
        end 
end
end
ylim([0 Yoff*(numGy+1)+200])
an=acg_metrics.acg_narrow(101:end,unit);
aw=acg_metrics.acg_wide(501:end,unit);
xLIN_narrow=linspace(1,400,101);
xLIN_wide = linspace(451,800,251);
plot(Xoff+xLIN_narrow,an/prctile(an,90)*200,'color','c','LineWidth',2.5);
plot([Xoff+xLIN_narrow(21) Xoff+xLIN_narrow(21)],[0 200],'c--','LineWidth',1)
plot([Xoff+xLIN_narrow(1) Xoff+xLIN_narrow(1)],[0 200],'c--','LineWidth',1)
plot(Xoff+xLIN_wide,aw(1:2:end)/prctile(aw,90)*200,'color','c','LineWidth',2.5);
plot([Xoff+xLIN_wide(51) Xoff+xLIN_wide(51)],[0 200],'c--','LineWidth',1)
plot([Xoff+xLIN_wide(1) Xoff+xLIN_wide(1)],[0 200],'c--','LineWidth',1)


plot([Xoff-3 Xoff-3],[-100 15*log2(FRForFeature(unit))]+100,'b','LineWidth',2)
plot([Xoff-3 Xoff],[100 100],'b--','LineWidth',1)

% try
%     xCenter = options.xCenterTransferFunction(round(mean(options.xCenter(options.currentComposition))));
%     xlim([(xCenter-5)*Xoff (xCenter+5)*Xoff])
% catch
% end
% end
% set(gca,'Position',[ 0.0895    0.0892    0.8112    0.8612])
set(gca,'XTickLabel','')
set(gca,'Position',[0 0 1 0.94])

axis off

title(sprintf('ses %04i unit# %04i',dayNumForEachUnit(unit),unit),'FontSize',22,'Color','k')
saveas(h,[shankStr '/videoFol/unit' num2str(unit,'%04i') '.png']);
% clf(h);
end
% close(h);
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