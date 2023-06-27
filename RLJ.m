function [posjadid] = RLJ(posc,posghabli,koljuv,shomare,Rprob,iff)
sumrr=zeros(1,size(posc,2));
for j=1:size(posc,2)
    for i=1: size(koljuv,1)
        sumrr(1,j)=sumrr(1,j)+koljuv(i).position(j);
    end  
end
meanrr=sumrr./size(koljuv,1);


 for i=1:size(posc,2)
    rand1=rand;
    if rand1<=Rprob
        alpha=posc-posghabli;
        beta=posghabli.*exp(-iff*rand*meanrr(i));
        if shomare<=50
        posjadid=posc+(beta-alpha);
        else
           posjadid=posc+(rand.*((rand.*beta)-alpha)); 
        end
    else
        posjadid=posc;
    end
 end
end
