function D2 = distfunXPPS(ZI, ZJ)
% ZJ(:,ZI==0)=0;
% ZI=repmat(ZI,size(ZJ,1),1);
% ZI(ZJ==0)=0;
[~,maxI]=max(abs(ZI),[],2);
[~,maxJ]=max(abs(ZJ),[],2);
% find common non-zero top ch ,( 0 blockage ) 
TempD = bsxfun(@max,abs(ZJ),abs(ZI));
TempN = abs(bsxfun(@minus,ZJ,ZI));
D1 = TempN./TempD;
D2=zeros(numel(maxJ),1);
for i = 1:numel(maxJ)
D2(i)=D1(i,maxJ(i))+D1(i,maxI);
end
end
%ZJ(:,ZI==0)=0;
% ZI=repmat(ZI,size(ZJ,1),1);
% ZI(ZJ==0)=0;
% [~,maxI]=max(abs(ZI),[],2);
% [~,maxJ]=max(abs(ZJ),[],2);
% % find common non-zero top ch ,( 0 blockage ) 
% TempD = max(abs(ZJ),abs(ZI));
% TempN = abs(ZJ-ZI);
% D1 = TempN./TempD;
% D2=zeros(numel(maxJ),1);
% for i = 1:numel(maxJ)
% D2(i)=D1(i,maxJ(i))+D1(i,maxI(i));
% end
