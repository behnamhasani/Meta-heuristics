function [posp1jadid,posp2jadid] = RLP(posp1,posp2,koljuv)
posc=posp1;
posp1jadid=zeros(1,size(posc,2));
posp2jadid=zeros(1,size(posc,2));
VarMin=-4;
VarMax=4;
sumrr=zeros(1,size(posc,2));
for j=1:size(posc,2)
    for i=1: size(koljuv,1)
        sumrr(1,j)=sumrr(1,j)+koljuv(i).position(j);
    end  
end
sumrr=sumrr+posp1+posp2;
meanrr=sumrr./(size(koljuv,1)+2);
for j=1:size(posc,2)
    
posp1jadid(1,j)=posp1(1,j)-(posp1(1,j)+exp(rand*(meanrr(j)-posp1(1,j))));

posp2jadid(1,j)=posp2(1,j)-(rand*(posp1(1,j)-exp(rand*(meanrr(j)-posp2(1,j)))));
if posp1jadid(1,j)<VarMin
            posp1jadid(1,j)=mod(posp1jadid(1,j),VarMin);
end
 if posp1jadid(1,j)>VarMax
             posp1jadid(1,j)=mod(posp1jadid(1,j),VarMax);
 end
 
 if posp2jadid(1,j)<VarMin
            posp2jadid(1,j)=mod(posp2jadid(1,j),VarMin);
end
 if posp2jadid(1,j)>VarMax
             posp2jadid(1,j)=mod(posp2jadid(1,j),VarMax);
 end


end
 


end