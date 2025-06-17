function D2 = distfunCosNaN(ZI, ZJ)
ZJ(:,isnan(ZI))=nan;
ZI=repmat(ZI,size(ZJ,1),1);
ZI(isnan(ZJ))=nan;

D2=zeros(size(ZJ,1),1);
for i = 1:size(ZJ,1)
    u=ZI(i,~isnan(ZI(i,:)));
    v=ZJ(i,~isnan(ZJ(i,:)));
D2(i)= (u*v')/norm(u)/norm(v);
end
D2=1-D2;
end