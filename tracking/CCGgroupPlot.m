function  CCGgroupPlot(HObj, event,composition,sesNumForEachUnit)
global ccg_PerSession withinSesValidCellId figCells
sesNumThisComp =sesNumForEachUnit(composition); 
N=histcounts(sesNumThisComp,0.5:1:max(sesNumForEachUnit)+0.5);
whichSesViolated = find(N>1);
for i = 1:numel(whichSesViolated)
figCells{i}=figure;
AcrossSesComp = composition(sesNumThisComp==whichSesViolated(i));
WithinThatSesComp = withinSesValidCellId(AcrossSesComp);
ccgSub = ccg_PerSession{whichSesViolated(i)}(:,WithinThatSesComp,WithinThatSesComp);
nU=numel(AcrossSesComp);
for j=1:nU-1
    for k = j:nU

     subplot(nU,nU,(j-1)*nU+k)
     ccg=squeeze(ccgSub(:,j,k));
     mid = round((numel(ccg)-1)/2);
     bar(ccg)
     left = sum(ccg(1:mid));
     right = sum(ccg(mid+1:end));
     Asym=max(left,right)/(min(left,right)+1); 
     title(['Asym = ' num2str(round(Asym,2))]);
     ylabel(num2str(AcrossSesComp(j)),'FontSize',12)
     xlabel(num2str(AcrossSesComp(k)),'FontSize',12)

            

    end
end
end
end
