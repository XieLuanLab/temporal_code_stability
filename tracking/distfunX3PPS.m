function D2 = distfunX3PPS(ZI, ZJ)
ZI=abs(ZI);
ZJ=abs(ZJ);

[~,sortI]=sort(ZI,2,'desc');
[~,sortJ]=sort(ZJ,2,'desc');

D2=zeros(size(ZJ,1),1);

for i = 1:size(ZJ,1)
  

jointTopCh = union(sortI(1:3) , sortJ(i,1:3));

% kI=polyfit(ZI(jointTopCh),ZJavg(i,jointTopCh),1);
% kJ=polyfit(ZJ(i,jointTopCh),ZJavg(i,jointTopCh),1);
% kI=Regression_fast([ones(sum(jointTopCh),1)
% ZI(jointTopCh)'],ZJavg(i,jointTopCh)'); too slow 
% kJ=Regression_fast([ones(sum(jointTopCh),1) ZJ(i,jointTopCh)'],ZJavg(i,jointTopCh)');


fI = ZI(jointTopCh);
fJ = ZJ(i,jointTopCh);

TempD = max(fI,fJ);
TempN = abs(fI-fJ);
D2(i)=sum(TempN./TempD)/numel(jointTopCh);
end

end


% %%
% ZI=abs(ZI);
% ZJ=abs(ZJ);
% ZJ(:,ZI==0)=0;
% ZI = repmat(ZI,size(ZJ,1),1);
% ZI(ZJ==0)=0;
% 
% [~,sortI]=sort(ZI,2,'desc');
% [~,sortJ]=sort(ZJ,2,'desc');
% 
% D2=zeros(size(ZJ,1),1);
% 
% for i = 1:size(ZJ,1)
%   
% 
% jointTopCh = union(sortI(i,1:3) , sortJ(i,1:3));
% 
% % kI=polyfit(ZI(jointTopCh),ZJavg(i,jointTopCh),1);
% % kJ=polyfit(ZJ(i,jointTopCh),ZJavg(i,jointTopCh),1);
% % kI=Regression_fast([ones(sum(jointTopCh),1)
% % ZI(jointTopCh)'],ZJavg(i,jointTopCh)'); too slow 
% % kJ=Regression_fast([ones(sum(jointTopCh),1) ZJ(i,jointTopCh)'],ZJavg(i,jointTopCh)');
% 
% 
% fI = ZI(i,jointTopCh);
% fJ = ZJ(i,jointTopCh);
% 
% TempD = max(fI,fJ);
% TempN = abs(fI-fJ);
% D2(i)=sum(TempN./TempD)/numel(jointTopCh);
% end